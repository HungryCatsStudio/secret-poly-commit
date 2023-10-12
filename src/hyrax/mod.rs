mod data_structures;
mod utils;
pub use data_structures::*;
mod tests;

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial};
use ark_std::{rand::RngCore, UniformRand};
use blake2::Blake2s256;
use core::marker::PhantomData;
use digest::Digest;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::linear_codes::utils::{
    inner_product, scalar_by_vector, vector_sum, IOPTranscript, Matrix,
};

use crate::{
    challenge::ChallengeGenerator,
    hyrax::utils::{flat_to_matrix_column_major, naive_chi, usize_to_bits},
    Error, LabeledCommitment, LabeledPolynomial, PolynomialCommitment,
};

/// String of bits used to seed the randomness during the setup function.
/// Note that the latter should never be used in production environments.
pub const PROTOCOL_NAME: &'static [u8] = b"Hyrax protocol";

/// Hyrax polynomial committment scheme:
/// A polynomial commitment scheme based on the hardness of the
/// discrete logarithm problem in prime-order groups. This is a
/// Fiat-Shamired version of the PCS described in the Hyrax paper
/// [[WTsTW17]][hyrax].
///
/// [hyrax]: https://eprint.iacr.org/2017/1132.pdf
///
/// * Modification note *
///
/// In the PCS contained in the cited article, the verifier never learns the
/// actual evaluation of the polynomial at the requested point, but is instead
/// convinced that a previously received Pedersen commitment is indeed a
/// commitment to said evaluation - this is what the SNARK proposed therein
/// necessitates. However, the Arkworks framework requies the verifier to
/// actually learn that value, which is why we have added the opening of
/// the commitment at the end of the protocol. This might not result in an
/// optimal non-hiding PCS, but we feel it is the most faithful adaptation of
/// original PCS that can be implemented with the current restrictions.
pub struct HyraxPC<
    // The curve used for Pedersen commitments (only EC groups are
    // supported as of now).
    G: AffineRepr,
    // TODO make generic or fix one type and remove this
    // S: CryptographicSponge,
> {
    _curve: PhantomData<G>,
}

// TODO so far everything is done with asserts instead of the Error
// types defined by the library. Is this okay?

// TODO use ark_std::cfg_iter! instead of iter() as it is now?

// TODO add "trusting" version of utils linear-algebra functions which do not check dimensions?

// TODO check if it makes sense to implement batch_check, batch_open, open_combinations or check_combinations

// TODO ********************************************************

// TODO document

impl<G: AffineRepr> HyraxPC<G> {
    fn pedersen_commit(
        key: &HyraxCommitterKey<G>,
        scalars: &[G::ScalarField],
        r: Option<G::ScalarField>,
        rng: Option<&mut dyn RngCore>,
    ) -> (G, G::ScalarField) {
        let r = r.unwrap_or(G::ScalarField::rand(rng.unwrap()));

        let mut scalars_ext = Vec::from(scalars);
        scalars_ext.push(r);

        let mut points_ext = key.com_key[0..scalars.len()].to_vec();
        points_ext.push(key.h);

        let scalars_bigint = scalars_ext
            .iter()
            .map(|s| s.into_bigint())
            .collect::<Vec<_>>();

        let com = <G::Group as VariableBaseMSM>::msm_bigint(&points_ext, &scalars_bigint);

        // TODO better way than with into? Difference AffineRep and idem::Group?
        (com.into(), r)
    }
}

type MLE<G: AffineRepr> = DenseMultilinearExtension<G::ScalarField>;

// TODO ********************************************************

impl<G: AffineRepr>
    PolynomialCommitment<
        G::ScalarField,
        DenseMultilinearExtension<G::ScalarField>,
        // Dummy sponge - required by the trait, not used in this implementation
        PoseidonSponge<G::ScalarField>,
    > for HyraxPC<G>
{
    type UniversalParams = HyraxUniversalParams<G>;
    type CommitterKey = HyraxCommitterKey<G>;
    type VerifierKey = HyraxVerifierKey<G>;
    type PreparedVerifierKey = HyraxPreparedVerifierKey<G>;
    type Commitment = HyraxCommitment<G>;
    type PreparedCommitment = HyraxPreparedCommitment<G>;
    type Randomness = HyraxRandomness<G::ScalarField>;
    type Proof = Vec<HyraxProof<G>>;
    type BatchProof = Vec<Self::Proof>;
    type Error = Error;

    /// Outputs mock universal parameters for the Hyrax polynomial commitment
    /// scheme. It does *not* return random keys across calls and should never
    /// be used in settings where security is required - it is only useful for
    /// testing. Furthermore, the point at infinity could possibly be part of
    /// the output, which sould not happen in an actual key.
    fn setup<R: RngCore>(
        max_degree: usize,
        num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        let n = num_vars.expect("Hyrax requires num_vars to be specified");

        assert_eq!(max_degree, 1, "Hyrax only supports multilinear polynomials");
        assert_eq!(
            n % 2,
            0,
            "Only polynomials with an even number of variables \
                    are supported in this implementation"
        );

        // Number of rows (or, equivalently, colums) of a square matrix
        // containing the coefficients of an n-variate ML polynomial
        let dim = 1 << n / 2;

        // The following block of code is largely taking from the IPA module
        // in this crate.
        let points: Vec<_> = ark_std::cfg_into_iter!(0u64..dim + 1)
            .map(|i| {
                let mut hash =
                    Blake2s256::digest([PROTOCOL_NAME, &i.to_le_bytes()].concat().as_slice());
                let mut p = G::from_random_bytes(&hash);
                let mut j = 0u64;
                while p.is_none() {
                    // PROTOCOL NAME, i, j
                    let mut bytes = b"Hyrax-protocol".to_vec();
                    bytes.extend(i.to_le_bytes());
                    bytes.extend(j.to_le_bytes());
                    hash = Blake2s256::digest(bytes.as_slice());
                    p = G::from_random_bytes(&hash);
                    j += 1;
                }
                let point = p.unwrap();
                point.mul_by_cofactor_to_group()
            })
            .collect();

        let mut points = G::Group::normalize_batch(&points);

        let h: G = points.pop().unwrap();

        Ok(HyraxUniversalParams {
            com_key: points,
            h,
            num_vars: n,
        })
    }

    fn trim(
        pp: &Self::UniversalParams,
        supported_degree: usize,
        supported_hiding_bound: usize,
        enforced_degree_bounds: Option<&[usize]>,
    ) -> Result<(Self::CommitterKey, Self::VerifierKey), Self::Error> {
        assert!(
            supported_degree == 1 && supported_hiding_bound == 1,
            "Hyrax only supports multilinear polynomials: \
            the passed degrees should be 1"
        );

        assert!(
            enforced_degree_bounds.is_none(),
            "Hyrax only supports multilinear polynomials: \
            enforced_degree_bounds should be `None`"
        );

        let HyraxUniversalParams {
            com_key,
            h,
            num_vars,
        } = pp.clone();

        let ck = HyraxCommitterKey {
            com_key,
            h,
            num_vars,
        };

        let vk: HyraxVerifierKey<G> = ck.clone();

        Ok((ck, vk))
    }

    /// Outputs a list of commitments to the passed polynomials
    fn commit<'a>(
        ck: &Self::CommitterKey,
        polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, MLE<G>>>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<
        (
            Vec<LabeledCommitment<Self::Commitment>>,
            Vec<Self::Randomness>,
        ),
        Self::Error,
    >
    where
        MLE<G>: 'a,
    {
        let mut coms = Vec::new();
        let mut rands = Vec::new();

        let rng_inner = rng.expect("Committing to polynomials requires a random generator");

        for l_poly in polynomials {
            let mut com_rands = Vec::new();

            let label = l_poly.label();
            let poly = l_poly.polynomial();

            assert_eq!(
                l_poly.degree_bound().unwrap_or(1),
                1,
                "Hyrax only supports ML polynomials: the degree bound should \
                be Some(1) or None"
            );

            assert!(
                l_poly.hiding_bound().is_none(),
                "Hiding bounds are not part of the Hyrax PCS"
            );

            let n = poly.num_vars();
            let dim = 1 << n / 2;

            assert!(
                n <= ck.num_vars,
                "Attempted to commit to a polynomial with {n} variables, but
                this key only supports up to {} variables",
                ck.num_vars
            );

            let m = flat_to_matrix_column_major(&poly.to_evaluations(), dim, dim);

            let row_coms = m
                .iter()
                .map(|row| {
                    let (c, r) = Self::pedersen_commit(ck, &row, None, Some(rng_inner));
                    com_rands.push(r);
                    c
                })
                .collect();

            let com = HyraxCommitment { row_coms };
            let l_comm = LabeledCommitment::new(label.to_string(), com, Some(1));

            coms.push(l_comm);
            rands.push(com_rands);
        }

        Ok((coms, rands))
    }

    fn open<'a>(
        ck: &Self::CommitterKey,
        labeled_polynomials: impl IntoIterator<Item = &'a LabeledPolynomial<G::ScalarField, MLE<G>>>,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a <DenseMultilinearExtension<G::ScalarField> as Polynomial<G::ScalarField>>::Point,
        // Not used and not generic on the cryptographic sponge S
        _opening_challenges: &mut ChallengeGenerator<
            G::ScalarField,
            PoseidonSponge<G::ScalarField>,
        >,
        rands: impl IntoIterator<Item = &'a Self::Randomness>,
        rng: Option<&mut dyn RngCore>,
    ) -> Result<Self::Proof, Self::Error>
    where
        Self::Commitment: 'a,
        Self::Randomness: 'a,
        MLE<G>: 'a,
    {
        // TODO is it safe to open several polynomials at once?
        // TODO is there a more efficient way to open several polynomials at once?
        //      can one e.g. share zs, ds..?

        let n = point.len();

        assert_eq!(
            n % 2,
            0,
            "Only points with an even number of variables \
            are supported in this implementation"
        );

        let dim = 1 << n / 2;

        let point_lower = &point[n / 2..];
        let point_upper = &point[..n / 2];

        // TODO this way to compute the bits is very inefficient
        let l: Vec<G::ScalarField> = (0..dim)
            .map(|idx| naive_chi(&usize_to_bits(idx, n / 2), point_lower))
            .collect();
        let r: Vec<G::ScalarField> = (0..dim)
            .map(|idx| naive_chi(&usize_to_bits(idx, n / 2), point_upper))
            .collect();

        let mut proofs = Vec::new();

        let rng_inner = rng.expect("Opening polynomials requires randomness");

        for (l_poly, (l_com, randomness)) in labeled_polynomials
            .into_iter()
            .zip(commitments.into_iter().zip(rands.into_iter()))
        {
            // TODO check if the poly was actually necessary
            let label = l_poly.label();
            assert_eq!(
                label,
                l_com.label(),
                "Mismatching labels: {label} and {}",
                l_com.label()
            );

            let poly = l_poly.polynomial();
            let com = l_com.commitment();

            // TODO chech num of vars matches n
            assert_eq!(
                poly.num_vars(),
                n,
                "The committed polynomial has {} variables, but \
                the point has {n} variables",
                poly.num_vars()
            );

            // Initialising the transcript
            let mut transcript: IOPTranscript<G::ScalarField> = IOPTranscript::new(b"transcript");

            // Absorbing public parameters
            transcript.append_serializable_element(b"public parameters", ck)?;

            // Absorbing the commitment to the polynomial
            transcript.append_serializable_element(b"commitment", &com.row_coms)?;

            // Absorbing the point
            transcript.append_serializable_element(b"point", point)?;

            // Commiting to the matrix formed by the polynomial coefficients
            let t_aux = flat_to_matrix_column_major(&poly.to_evaluations(), dim, dim);
            let t = Matrix::new_from_rows(t_aux);

            let lt = t.row_mul(&l);

            let r_lt = l
                .iter()
                .zip(randomness.iter())
                .map(|(l, r)| *l * r)
                .sum::<G::ScalarField>();

            let eval = inner_product(&lt, &r);

            // Singleton commit
            let (com_eval, r_eval) = Self::pedersen_commit(ck, &[eval], None, Some(rng_inner));

            // ******** Dot product argument ********
            // Appendix A.2 in the reference article

            let d: Vec<G::ScalarField> =
                (0..dim).map(|_| G::ScalarField::rand(rng_inner)).collect();

            let b = inner_product(&r, &d);

            // Multi-commit
            let (com_d, r_d) = Self::pedersen_commit(ck, &d, None, Some(rng_inner));

            // Singleton commit
            let (com_b, r_b) = Self::pedersen_commit(ck, &[b], None, Some(rng_inner));

            // Absorbing the commitment to the evaluation
            transcript.append_serializable_element(b"com_eval", &com_eval)?;

            // Absorbing the two auxiliary commitments
            transcript.append_serializable_element(b"com_d", &com_d)?;
            transcript.append_serializable_element(b"com_b", &com_b)?;

            // Receive the random challenge c from the verifier, i.e. squeeze
            // it from the transcript.
            let c = transcript.get_and_append_challenge(b"c").unwrap();

            let z = vector_sum(&d, &scalar_by_vector(c, &lt));
            let z_d = c * r_lt + r_d;
            let z_b = c * r_eval + r_b;

            // ******** Opening ********
            // This is *not* part of the Hyrax PCS as described in the reference
            // article. Cf. the "Modification note" at the beginning of this file.
            // From the prover's perspective, opening amounts to adding r_eval to
            // the proof.

            proofs.push(HyraxProof {
                com_eval,
                com_d,
                com_b,
                z,
                z_d,
                z_b,
                r_eval,
            });
        }

        Ok(proofs)
    }

    fn check<'a>(
        vk: &Self::VerifierKey,
        commitments: impl IntoIterator<Item = &'a LabeledCommitment<Self::Commitment>>,
        point: &'a <DenseMultilinearExtension<G::ScalarField> as Polynomial<G::ScalarField>>::Point,
        values: impl IntoIterator<Item = G::ScalarField>,
        proof: &Self::Proof,
        // Not used and not generic on the cryptographic sponge S
        _opening_challenges: &mut ChallengeGenerator<
            G::ScalarField,
            PoseidonSponge<G::ScalarField>,
        >,
        _rng: Option<&mut dyn RngCore>,
    ) -> Result<bool, Self::Error>
    where
        Self::Commitment: 'a,
    {
        let n = point.len();

        assert_eq!(
            n % 2,
            0,
            "Only points with an even number of variables \
            are supported in this implementation"
        );

        let dim = 1 << n / 2;

        let point_lower = &point[n / 2..];
        let point_upper = &point[..n / 2];

        // TODO this way to compute the bits is very inefficient
        let l: Vec<G::ScalarField> = (0..dim)
            .map(|idx| naive_chi(&usize_to_bits(idx, n / 2), point_lower))
            .collect();
        let r: Vec<G::ScalarField> = (0..dim)
            .map(|idx| naive_chi(&usize_to_bits(idx, n / 2), point_upper))
            .collect();

        for (com, (claim, h_proof)) in commitments
            .into_iter()
            .zip(values.into_iter().zip(proof.into_iter()))
        {
            let row_coms = &com.commitment().row_coms;

            // extract each field from h_proof
            let HyraxProof {
                com_eval,
                com_d,
                com_b,
                z,
                z_d,
                z_b,
                r_eval,
            } = h_proof;

            // TODO chech num of vars matches n
            assert_eq!(
                row_coms.len(),
                1 << n / 2,
                "The commitment should have 2^(n/2) = has {} entries, but \
                it has {} instead",
                1 << n / 2,
                row_coms.len()
            );

            // TODO change to multi-exponentiation OR directly compute as the commitment to LT?
            let t_prime: G = row_coms
                .iter()
                .zip(l.iter())
                .map(|(e, s)| e.mul(s))
                .sum::<G::Group>()
                .into();

            // Construct transcript and squeeze the challenge c from it

            let mut transcript: IOPTranscript<G::ScalarField> = IOPTranscript::new(b"transcript");

            // Absorbing public parameters
            transcript.append_serializable_element(b"public parameters", vk)?;

            // Absorbing the commitment to the polynomial
            transcript.append_serializable_element(b"commitment", row_coms)?;

            // Absorbing the point
            transcript.append_serializable_element(b"point", point)?;

            // Absorbing the commitment to the evaluation
            transcript.append_serializable_element(b"com_eval", com_eval)?;

            // Absorbing the two auxiliary commitments
            transcript.append_serializable_element(b"com_d", com_d)?;
            transcript.append_serializable_element(b"com_b", com_b)?;

            // Receive the random challenge c from the verifier, i.e. squeeze
            // it from the transcript.
            let c = transcript.get_and_append_challenge(b"c").unwrap();

            // First check
            let com_z_zd = Self::pedersen_commit(vk, &z, Some(*z_d), None).0;
            if com_z_zd != (t_prime.mul(c) + com_d).into() {
                return Ok(false);
            }

            // Second check
            let com_dp = Self::pedersen_commit(vk, &[inner_product(&r, &z)], Some(*z_b), None).0;
            // TODO clarify why into() is needed
            if com_dp != (com_eval.mul(c) + com_b).into() {
                return Ok(false);
            }

            // Third check: opening
            let exp = Self::pedersen_commit(vk, &[claim], Some(*r_eval), None).0;
            if *com_eval != exp {
                return Ok(false);
            }
        }

        Ok(true)
    }
}