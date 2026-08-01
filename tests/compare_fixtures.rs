#[test]
fn compare_matches_expected_fixture_outputs() {
    hap_rs::verification::verify_rust_outputs().expect("rust fixture verification failed");
}
