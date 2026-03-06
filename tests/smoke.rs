use rotamer::*;

mod common;
use common::*;

macro_rules! smoke {
    ($T:ident, 0, 0, $n:literal) => {
        paste::paste! {
            #[test]
            fn [<smoke_ $T:lower>]() {
                let coords = $T::build(N, CA, C);
                let s = coords.as_slice();
                assert_eq!(
                    s.len(),
                    $n,
                    "{}: expected {} atoms, got {}",
                    stringify!($T),
                    $n,
                    s.len()
                );
                assert_coords_valid(s, stringify!($T));
            }
        }
    };
    ($T:ident, $nc:literal, 0, $n:literal) => {
        paste::paste! {
            #[test]
            fn [<smoke_ $T:lower>]() {
                let coords = $T::build(N, CA, C, [0.0_f32; $nc]);
                let s = coords.as_slice();
                assert_eq!(
                    s.len(),
                    $n,
                    "{}: expected {} atoms, got {}",
                    stringify!($T),
                    $n,
                    s.len()
                );
                assert_coords_valid(s, stringify!($T));
            }
        }
    };
    ($T:ident, $nc:literal, $np:literal, $n:literal) => {
        paste::paste! {
            #[test]
            fn [<smoke_ $T:lower>]() {
                let coords = $T::build(N, CA, C, [0.0_f32; $nc], [0.0_f32; $np]);
                let s = coords.as_slice();
                assert_eq!(
                    s.len(),
                    $n,
                    "{}: expected {} atoms, got {}",
                    stringify!($T),
                    $n,
                    s.len()
                );
                assert_coords_valid(s, stringify!($T));
            }
        }
    };
}

for_all_sidechains!(smoke);
