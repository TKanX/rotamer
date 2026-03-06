use rotamer::*;

macro_rules! check_constants {
    ($T:ident, $nc:literal, $np:literal, $n:literal) => {
        paste::paste! {
            #[test]
            fn [<constants_ $T:lower>]() {
                assert_eq!(
                    <$T as Sidechain>::N_CHI,
                    $nc,
                    "{}: N_CHI mismatch",
                    stringify!($T)
                );
                assert_eq!(
                    <$T as Sidechain>::N_POLAR_H,
                    $np,
                    "{}: N_POLAR_H mismatch",
                    stringify!($T)
                );
                assert_eq!(
                    <<$T as Sidechain>::Coords as SidechainCoords>::N,
                    $n,
                    "{}: N mismatch",
                    stringify!($T)
                );
                let name = <$T as Sidechain>::NAME;
                assert_eq!(name.len(), 3, "{}: NAME '{}' is not 3 chars", stringify!($T), name);
                assert!(
                    name.bytes().all(|b| b.is_ascii_uppercase()),
                    "{}: NAME '{}' should be all-uppercase ASCII",
                    stringify!($T),
                    name
                );
                assert_eq!(
                    name,
                    stringify!($T).to_uppercase(),
                    "{}: NAME '{}' does not match type name",
                    stringify!($T),
                    name
                );
            }
        }
    };
}

for_all_sidechains!(check_constants);
