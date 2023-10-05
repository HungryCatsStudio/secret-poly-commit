#[cfg(test)]
mod tests {

    use crate::ligero::LinearCodePCS;
    use crate::{
        challenge::ChallengeGenerator,
        ligero::{utils::*, LigeroPCUniversalParams, MultilinearLigero, PolynomialCommitment},
        LabeledPolynomial,
    };
    use ark_bls12_377::Fq;
    use ark_bls12_377::Fr;
    use ark_bls12_381::Fr as Fr381;
    use ark_crypto_primitives::sponge::poseidon::PoseidonConfig;
    use ark_crypto_primitives::sponge::CryptographicSponge;
    use ark_crypto_primitives::{
        crh::{pedersen, sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
        merkle_tree::{ByteDigestConverter, Config},
        sponge::poseidon::PoseidonSponge,
    };
    use ark_ff::{Field, PrimeField};
    use ark_poly::evaluations::multivariate::{MultilinearExtension, SparseMultilinearExtension};
    use ark_std::test_rng;
    use blake2::Blake2s256;
    use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

    #[derive(Clone)]
    pub(super) struct Window4x256;
    impl pedersen::Window for Window4x256 {
        const WINDOW_SIZE: usize = 4;
        const NUM_WINDOWS: usize = 256;
    }

    type LeafH = Sha256;
    type CompressH = Sha256;

    struct MerkleTreeParams;

    impl Config for MerkleTreeParams {
        type Leaf = [u8];

        type LeafDigest = <LeafH as CRHScheme>::Output;
        type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
        type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

        type LeafHash = LeafH;
        type TwoToOneHash = CompressH;
    }

    type MTConfig = MerkleTreeParams;
    type Sponge = PoseidonSponge<Fr>;

    type LigeroPCS = LinearCodePCS<
        MultilinearLigero<Fr, MTConfig, Blake2s256, Sponge, SparseMultilinearExtension<Fr>>,
        Fr,
        SparseMultilinearExtension<Fr>,
        Sponge,
        MTConfig,
        Blake2s256,
    >;
    type LigeroPcsF<F> = LinearCodePCS<
        MultilinearLigero<F, MTConfig, Blake2s256, Sponge, SparseMultilinearExtension<F>>,
        F,
        SparseMultilinearExtension<F>,
        Sponge,
        MTConfig,
        Blake2s256,
    >;

    fn rand_poly<Fr: PrimeField>(
        _: usize,
        num_vars: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> SparseMultilinearExtension<Fr> {
        match num_vars {
            Some(n) => SparseMultilinearExtension::rand(n, rng),
            None => unimplemented!(), // should not happen in ML case!
        }
    }

    fn constant_poly<Fr: PrimeField>(
        _: usize,
        num_vars: Option<usize>,
        rng: &mut ChaCha20Rng,
    ) -> SparseMultilinearExtension<Fr> {
        // f1 = (1-x1)(1-x2)(1-x3)(1-x5)[(1-x6)*x4 + 2(1-x4)*x6]
        match num_vars {
            Some(n) => {
                let points = vec![(1, Fr::rand(rng))];
                SparseMultilinearExtension::from_evaluations(n, &points)
            }
            None => unimplemented!(), // should not happen in ML case!
        }
    }

    // TODO: replace by https://github.com/arkworks-rs/crypto-primitives/issues/112.
    fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
        let full_rounds = 8;
        let partial_rounds = 31;
        let alpha = 17;

        let mds = vec![
            vec![F::one(), F::zero(), F::one()],
            vec![F::one(), F::one(), F::zero()],
            vec![F::zero(), F::one(), F::one()],
        ];

        let mut v = Vec::new();
        let mut ark_rng = test_rng();

        for _ in 0..(full_rounds + partial_rounds) {
            let mut res = Vec::new();

            for _ in 0..3 {
                res.push(F::rand(&mut ark_rng));
            }
            v.push(res);
        }
        let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
        PoseidonSponge::new(&config)
    }

    #[test]
    fn test_construction() {
        let mut rng = &mut test_rng();
        // just to make sure we have the right degree given the FFT domain for our field
        let leaf_hash_params = <LeafH as CRHScheme>::setup(&mut rng).unwrap();
        let two_to_one_params = <CompressH as TwoToOneCRHScheme>::setup(&mut rng)
            .unwrap()
            .clone();
        let check_well_formedness = true;

        let pp: LigeroPCUniversalParams<Fr, MTConfig> = LigeroPCUniversalParams::new(
            128,
            4,
            check_well_formedness,
            leaf_hash_params,
            two_to_one_params,
        );

        let (ck, vk) = LigeroPCS::trim(&pp, 0, 0, None).unwrap();

        let rand_chacha = &mut ChaCha20Rng::from_rng(test_rng()).unwrap();
        let labeled_poly = LabeledPolynomial::new(
            "test".to_string(),
            rand_poly(1, Some(5), rand_chacha),
            Some(5),
            Some(5),
        );

        let mut test_sponge = test_sponge::<Fr>();
        let (c, rands) = LigeroPCS::commit(&ck, &[labeled_poly.clone()], None).unwrap();

        let point = rand_point(Some(5), rand_chacha);

        let value = labeled_poly.evaluate(&point);

        let mut challenge_generator: ChallengeGenerator<Fr, PoseidonSponge<Fr>> =
            ChallengeGenerator::new_univariate(&mut test_sponge);

        let proof = LigeroPCS::open(
            &ck,
            &[labeled_poly],
            &c,
            &point,
            &mut (challenge_generator.clone()),
            &rands,
            None,
        )
        .unwrap();
        assert!(LigeroPCS::check(
            &vk,
            &c,
            &point,
            [value],
            &proof,
            &mut challenge_generator,
            None
        )
        .unwrap());
    }

    #[test]
    fn test_calculate_t_with_good_parameters() {
        assert!(calculate_t::<Fq>(128, 4, 2_usize.pow(32)).unwrap() < 200);
        assert!(calculate_t::<Fq>(256, 4, 2_usize.pow(32)).unwrap() < 400);
    }

    #[test]
    fn test_calculate_t_with_bad_parameters() {
        calculate_t::<Fq>((Fq::MODULUS_BIT_SIZE - 60) as usize, 4, 2_usize.pow(60)).unwrap_err();
        calculate_t::<Fq>(400, 4, 2_usize.pow(32)).unwrap_err();
    }

    fn rand_point<F: Field>(num_vars: Option<usize>, rng: &mut ChaCha20Rng) -> Vec<F> {
        match num_vars {
            Some(n) => (0..n).map(|_| F::rand(rng)).collect(),
            None => unimplemented!(), // should not happen!
        }
    }

    #[test]
    fn single_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, LigeroPcsF<Fr381>, _>(
            Some(10),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn constant_poly_test() {
        use crate::tests::*;
        single_poly_test::<_, _, LigeroPCS, _>(
            Some(10),
            constant_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        single_poly_test::<_, _, LigeroPcsF<Fr381>, _>(
            Some(5),
            constant_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
    }

    #[test]
    fn full_end_to_end_test() {
        use crate::tests::*;
        full_end_to_end_test::<_, _, LigeroPCS, _>(
            Some(8),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_test::<_, _, LigeroPcsF<Fr381>, _>(
            Some(3),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn single_equation_test() {
        use crate::tests::*;
        single_equation_test::<_, _, LigeroPCS, _>(
            Some(10),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        single_equation_test::<_, _, LigeroPcsF<Fr381>, _>(
            Some(5),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn two_equation_test() {
        use crate::tests::*;
        two_equation_test::<_, _, LigeroPCS, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        two_equation_test::<_, _, LigeroPcsF<Fr381>, _>(
            Some(10),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }

    #[test]
    fn full_end_to_end_equation_test() {
        use crate::tests::*;
        full_end_to_end_equation_test::<_, _, LigeroPCS, _>(
            Some(5),
            rand_poly::<Fr>,
            rand_point::<Fr>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-377");
        println!("Finished bls12-377");
        full_end_to_end_equation_test::<_, _, LigeroPcsF<Fr381>, _>(
            Some(8),
            rand_poly::<Fr381>,
            rand_point::<Fr381>,
            poseidon_sponge_for_test,
        )
        .expect("test failed for bls12-381");
        println!("Finished bls12-381");
    }
}
