use super::utils::reed_solomon;
use super::{LigeroPCParams, LinCodeInfo, LinearEncode};

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{merkle_tree::Config, sponge::CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;
use ark_std::marker::PhantomData;
use ark_std::vec::Vec;
use digest::Digest;

mod tests;

/// The univariate Ligero polynomial commitment scheme based on [[Ligero]][ligero].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [ligero]: https://eprint.iacr.org/2022/1608.pdf
pub struct UnivariateLigero<
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
    H: CRHScheme,
> {
    _phantom: PhantomData<(F, C, D, S, P, H)>,
}

impl<F, C, D, S, P, H> LinearEncode<F, C, D, P, H> for UnivariateLigero<F, C, D, S, P, H>
where
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: DenseUVPolynomial<F>,
    P::Point: Into<F>,
    H: CRHScheme,
{
    type LinCodePCParams = LigeroPCParams<F, C, H>;

    fn setup<R>(
        _rng: &mut R,
        leaf_hash_params: <<C as Config>::LeafHash as CRHScheme>::Parameters,
        two_to_one_params: <<C as Config>::TwoToOneHash as TwoToOneCRHScheme>::Parameters,
        col_hash_params: H::Parameters,
    ) -> Self::LinCodePCParams {
        Self::LinCodePCParams::new(
            128,
            4,
            true,
            leaf_hash_params,
            two_to_one_params,
            col_hash_params,
        )
    }

    fn encode(msg: &[F], param: &Self::LinCodePCParams) -> Vec<F> {
        reed_solomon(msg, param.rho_inv().0)
    }

    /// For a univariate polynomial, we simply return the list of coefficients.ś
    fn poly_repr(polynomial: &P) -> Vec<F> {
        polynomial.coeffs().to_vec()
    }

    fn point_to_vec(point: P::Point) -> Vec<F> {
        vec![point]
    }

    /// Compute out = [1, z, z^2, ..., z^(n_cols_1)]
    fn tensor(z: &F, left: usize, right: usize) -> (Vec<F>, Vec<F>) {
        let mut left_out = Vec::with_capacity(left);
        let mut pow_a = F::one();
        for _ in 0..left {
            left_out.push(pow_a);
            pow_a *= z;
        }

        let mut right_out = Vec::with_capacity(right);
        let mut pow_b = F::one();
        for _ in 0..right {
            right_out.push(pow_b);
            pow_b *= pow_a;
        }

        (left_out, right_out)
    }
}
