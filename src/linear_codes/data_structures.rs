use super::utils::SprsMat;
use crate::{
    PCCommitment, PCPreparedCommitment, PCPreparedVerifierKey, PCRandomness, PCVerifierKey,
};
use ark_crypto_primitives::{
    crh::CRHScheme,
    merkle_tree::{Config, LeafParam, Path, TwoToOneParam},
};
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::marker::PhantomData;
use ark_std::rand::RngCore;
use ark_std::vec::Vec;

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
/// The public parameters for Ligero PCS.
pub struct LigeroPCParams<F: PrimeField, C: Config, H: CRHScheme> {
    pub(crate) _field: PhantomData<F>,
    /// The security parameter
    pub(crate) sec_param: usize,
    /// The inverse of the code rate.
    pub(crate) rho_inv: usize,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub(crate) leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub(crate) two_to_one_params: TwoToOneParam<C>,
    #[derivative(Debug = "ignore")]
    pub(crate) col_hash_params: H::Parameters,
}

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone(bound = ""), Debug(bound = ""))]
/// The public parameters for Ligero PCS.
pub struct BrakedownPCParams<F: PrimeField, C: Config, H: CRHScheme> {
    /// The security parameter
    pub(crate) sec_param: usize,
    /// alpha in the paper
    pub(crate) alpha: (usize, usize),
    /// beta in the paper
    pub(crate) beta: (usize, usize),
    /// The inverse of the code rate.
    pub(crate) rho_inv: (usize, usize),
    /// Size of the base case to encode with RS
    pub(crate) base_len: usize,
    /// Number of rows in polynomial matrix
    pub(crate) n: usize,
    /// Number of columns in polynomial matrix
    pub(crate) m: usize,
    /// Length of codeword
    pub(crate) m_ext: usize,
    /// Size and sparsity of all of matrices A
    pub(crate) a_dims: Vec<(usize, usize, usize)>,
    /// Size and sparsity of all of matrices B
    pub(crate) b_dims: Vec<(usize, usize, usize)>,
    /// Start indices of the middle chunks during encoding
    pub(crate) start: Vec<usize>,
    /// End indices of the middle chunks during encoding
    pub(crate) end: Vec<usize>,
    /// Matrices
    pub(crate) a_mats: Vec<SprsMat<F>>,
    /// Matrices
    pub(crate) b_mats: Vec<SprsMat<F>>,
    /// This is a flag which determines if the random linear combination is done.
    pub(crate) check_well_formedness: bool,
    /// Parameters for hash function of Merkle tree leaves
    #[derivative(Debug = "ignore")]
    pub(crate) leaf_hash_params: LeafParam<C>,
    /// Parameters for hash function of Merke tree combining two nodes into one
    #[derivative(Debug = "ignore")]
    pub(crate) two_to_one_params: TwoToOneParam<C>,
    #[derivative(Debug = "ignore")]
    pub(crate) col_hash_params: H::Parameters,
}

pub(crate) type LinCodePCPreparedVerifierKey = ();

impl<Unprepared: PCVerifierKey> PCPreparedVerifierKey<Unprepared> for LinCodePCPreparedVerifierKey {
    fn prepare(_vk: &Unprepared) -> Self {}
}
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub(crate) struct Metadata {
    pub(crate) n_rows: usize,
    pub(crate) n_cols: usize,
    pub(crate) n_ext_cols: usize,
}

/// The commitment to a polynomial is a root of the merkle tree,
/// where each node is a hash of the column of the encoded coefficient matrix U.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCCommitment<C: Config> {
    // number of rows resp. columns of the square matrix containing the coefficients of the polynomial
    pub(crate) metadata: Metadata,
    pub(crate) root: C::InnerDigest,
}

impl<C: Config> PCCommitment for LinCodePCCommitment<C> {
    fn empty() -> Self {
        LinCodePCCommitment::default()
    }

    fn has_degree_bound(&self) -> bool {
        false
    }
}

pub(crate) type LinCodePCPreparedCommitment<C> = LinCodePCCommitment<C>;

impl<Unprepared: PCCommitment, C: Config> PCPreparedCommitment<Unprepared>
    for LinCodePCPreparedCommitment<C>
{
    fn prepare(_cm: &Unprepared) -> Self {
        LinCodePCPreparedCommitment::default()
    }
}

pub(crate) type LinCodePCRandomness = ();

impl PCRandomness for LinCodePCRandomness {
    fn empty() -> Self {
        unimplemented!()
    }

    fn rand<R: RngCore>(
        _num_queries: usize,
        _has_degree_bound: bool,
        _num_vars: Option<usize>,
        _rng: &mut R,
    ) -> Self {
        unimplemented!()
    }
}

/// Proof of an individual linear code well-formedness check or opening
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub(crate) struct LinCodePCProofSingle<F, C>
where
    F: PrimeField,
    C: Config,
{
    /// For each of the indices in q, `paths` contains the path from the root of the merkle tree to the leaf
    pub(crate) paths: Vec<Path<C>>,

    /// v, s.t. E(v) = w
    pub(crate) v: Vec<F>,

    pub(crate) columns: Vec<Vec<F>>,
}

/// The Proof type for linear code PCS, which amounts to an array of individual proofs
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Default(bound = ""), Clone(bound = ""), Debug(bound = ""))]
pub struct LinCodePCProof<F, C>
where
    F: PrimeField,
    C: Config,
{
    pub(crate) opening: LinCodePCProofSingle<F, C>,
    pub(crate) well_formedness: Option<Vec<F>>,
}

// Multiple poly at one point
pub(crate) type LPCPArray<F, C> = Vec<LinCodePCProof<F, C>>;
