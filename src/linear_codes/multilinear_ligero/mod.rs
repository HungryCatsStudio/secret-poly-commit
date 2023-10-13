use core::usize;

use ark_crypto_primitives::crh::{CRHScheme, TwoToOneCRHScheme};
use ark_crypto_primitives::{merkle_tree::Config, sponge::CryptographicSponge};
use ark_ff::{FftField, PrimeField};
use ark_poly::{MultilinearExtension, Polynomial};
use ark_std::borrow::Borrow;
use ark_std::log2;
use ark_std::marker::PhantomData;
use ark_std::vec::Vec;

use digest::Digest;

use super::utils::reed_solomon;
use super::{LigeroPCParams, LinCodeInfo, LinearEncode};

mod tests;

/// The multilinear Ligero polynomial commitment scheme based on [[Ligero]][ligero].
/// The scheme defaults to the naive batching strategy.
///
/// Note: The scheme currently does not support hiding.
///
/// [ligero]: https://eprint.iacr.org/2022/1608.pdf
pub struct MultilinearLigero<
    F: PrimeField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: MultilinearExtension<F>,
    H: CRHScheme,
> {
    _phantom: PhantomData<(F, C, D, S, P, H)>,
}

impl<F, C, D, S, P, H> LinearEncode<F, P, C, D, H> for MultilinearLigero<F, C, D, S, P, H>
where
    F: PrimeField + FftField,
    C: Config,
    D: Digest,
    S: CryptographicSponge,
    P: MultilinearExtension<F>,
    Vec<u8>: Borrow<C::Leaf>,
    <P as Polynomial<F>>::Point: Into<Vec<F>>,
    H: CRHScheme,
{
    type LinCodePCParams = LigeroPCParams<F, C, H>;

    fn setup(
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

    fn poly_repr(polynomial: &P) -> Vec<F> {
        polynomial.to_evaluations()
    }

    fn point_to_vec(point: <P as Polynomial<F>>::Point) -> Vec<F> {
        point.into()
    }

    fn tensor(
        point: &<P as Polynomial<F>>::Point,
        left_len: usize,
        _right_len: usize,
    ) -> (Vec<F>, Vec<F>) {
        let point: Vec<F> = Self::point_to_vec(point.clone());

        let split = log2(left_len) as usize;
        let left = &point[..split];
        let right = &point[split..];
        (tensor_inner(left), tensor_inner(right))
    }
}

fn tensor_inner<F: PrimeField>(values: &[F]) -> Vec<F> {
    let one = F::one();
    let anti_values: Vec<F> = values.iter().map(|v| one - *v).collect();

    let mut layer: Vec<F> = vec![one];

    for i in 0..values.len() {
        let mut new_layer = Vec::new();
        for v in &layer {
            new_layer.push(*v * anti_values[i]);
        }
        for v in &layer {
            new_layer.push(*v * values[i]);
        }
        layer = new_layer;
    }

    layer
}
