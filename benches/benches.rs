#![cfg(feature = "benches")]
use ark_ec::AffineRepr;
use ark_poly::DenseMultilinearExtension;
use blake2::Blake2s256;
use criterion::{criterion_group, criterion_main, Criterion};

use ark_crypto_primitives::{
    crh::{sha256::Sha256, CRHScheme, TwoToOneCRHScheme},
    merkle_tree::{ByteDigestConverter, Config},
    sponge::poseidon::PoseidonSponge,
};

use ark_poly_commit::{
    bench_templates::{bench_pcs_method, commit, open, verify, MLE},
    hyrax::HyraxPC,
    linear_codes::{
        FieldToBytesColHasher, LeafIdentityHasher, LinearCodePCS, MultilinearBrakedown,
        MultilinearLigero,
    },
};

use ark_bls12_381::{Fr as Fr381, G1Affine as G1Affine381};
use ark_bn254::{Fr as Fr254, G1Affine as G1Affine254};

// Hyrax type alias
type Hyrax<G> = HyraxPC<G, DenseMultilinearExtension<<G as AffineRepr>::ScalarField>>;

struct MerkleTreeParams;
type LeafH = LeafIdentityHasher;
type CompressH = Sha256;
impl Config for MerkleTreeParams {
    type Leaf = Vec<u8>;

    type LeafDigest = <LeafH as CRHScheme>::Output;
    type LeafInnerDigestConverter = ByteDigestConverter<Self::LeafDigest>;
    type InnerDigest = <CompressH as TwoToOneCRHScheme>::Output;

    type LeafHash = LeafH;
    type TwoToOneHash = CompressH;
}

type MTConfig = MerkleTreeParams;
type Sponge<F> = PoseidonSponge<F>;
type ColHasher<F> = FieldToBytesColHasher<F, Blake2s256>;

// Ligero type alias
type Ligero<F> = LinearCodePCS<
    MultilinearLigero<F, MTConfig, Sponge<F>, MLE<F>, ColHasher<F>>,
    F,
    MLE<F>,
    Sponge<F>,
    MTConfig,
    ColHasher<F>,
>;

// Brakedown type alias
type Brakedown<F> = LinearCodePCS<
    MultilinearBrakedown<F, MTConfig, Sponge<F>, MLE<F>, ColHasher<F>>,
    F,
    MLE<F>,
    Sponge<F>,
    MTConfig,
    ColHasher<F>,
>;

const MIN_NUM_VARS: usize = 12;
const MAX_NUM_VARS: usize = 25;

/*************** Instantiating target functions ***************/

fn hyrax_bn254(c: &mut Criterion) {
    bench_pcs_method::<_, Hyrax<G1Affine254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(4).collect(),
        "commit_hyrax_range_BN_254",
        commit::<_, Hyrax<G1Affine254>>,
    );
    bench_pcs_method::<_, Hyrax<G1Affine254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(4).collect(),
        "open_hyrax_range_BN_254",
        open::<_, Hyrax<G1Affine254>>,
    );

    bench_pcs_method::<_, Hyrax<G1Affine254>>(
        c,
        (MIN_NUM_VARS..MAX_NUM_VARS).step_by(4).collect(),
        "verify_hyrax_range_BN_254",
        verify::<_, Hyrax<G1Affine254>>,
    );
}

criterion_group! {
    name = hyrax_benches;
    config = Criterion::default();
    targets = hyrax_bn254
}

criterion_main!(hyrax_benches);
