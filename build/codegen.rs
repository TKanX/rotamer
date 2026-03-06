use std::f64::consts::PI;
use std::fs;
use std::path::PathBuf;

use proc_macro2::{Ident, Literal, Span, TokenStream};
use quote::{format_ident, quote};

use crate::residues::{ALL_RESIDUES, AtomSpec, ResidueSpec, TorsionSrc};

fn field_ident(pdb: &str) -> Ident {
    Ident::new(&pdb.to_ascii_lowercase(), Span::call_site())
}

fn f32_tok(v: f64) -> TokenStream {
    format!("{}_f32", v as f32)
        .parse::<TokenStream>()
        .expect("f32 literal")
}

fn trig_pair_tok(deg: f64) -> TokenStream {
    let rad = deg * PI / 180.0;
    let cos = f32_tok(rad.cos());
    let sin = f32_tok(rad.sin());
    quote! { TrigPair { cos: #cos, sin: #sin } }
}

pub fn generate() {
    let out_dir = PathBuf::from(std::env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"));
    let dest = out_dir.join("generated.rs");

    let blocks: Vec<TokenStream> = ALL_RESIDUES.iter().map(emit_residue).collect();
    let for_all = emit_for_all_sidechains_macro();
    let tokens = quote! { #(#blocks)* #for_all };

    fs::write(
        &dest,
        prettyplease::unparse(&syn::parse2(tokens).expect("invalid Rust")),
    )
    .expect("failed to write generated.rs");
}

fn emit_residue(res: &ResidueSpec) -> TokenStream {
    validate_refs(res);

    let ty = Ident::new(res.type_name, Span::call_site());
    let coords = format_ident!("{}Coords", res.type_name);
    let name = res.name;
    let n = res.atoms.len();
    let n_chi = res.n_chi;
    let n_ph = res.n_polar_h;
    let fields: Vec<Ident> = res.atoms.iter().map(|a| field_ident(a.name)).collect();

    let assert_msg = format!("{coords} layout must match [Vec3; {n}]");
    let struct_block = if n == 0 {
        quote! {
            #[repr(C)]
            #[derive(Debug, Clone, Copy, PartialEq)]
            pub struct #coords {}
        }
    } else {
        quote! {
            #[repr(C)]
            #[derive(Debug, Clone, Copy, PartialEq)]
            pub struct #coords {
                #(pub #fields: Vec3,)*
            }
            const _: () = assert!(
                core::mem::size_of::<#coords>() == #n * core::mem::size_of::<Vec3>(),
                #assert_msg
            );
        }
    };

    let as_slice_body = if n == 0 {
        quote! { &[] }
    } else {
        quote! {
            #[allow(clippy::undocumented_unsafe_blocks)]
            unsafe {
                core::slice::from_raw_parts(self as *const Self as *const Vec3, #n)
            }
        }
    };

    let coords_impl = quote! {
        impl SidechainCoords for #coords {
            const N: usize = #n;
            #[inline(always)]
            fn as_slice(&self) -> &[Vec3] {
                #as_slice_body
            }
        }
    };

    let sidechain_impl = quote! {
        impl Sidechain for #ty {
            const N_CHI: usize = #n_chi;
            const N_POLAR_H: usize = #n_ph;
            const NAME: &'static str = #name;
            type Coords = #coords;
        }
    };

    let build_impl = emit_build(res, &coords, &ty);

    quote! {
        #struct_block
        #coords_impl
        #sidechain_impl
        #build_impl
    }
}

fn emit_build(res: &ResidueSpec, coords: &Ident, ty: &Ident) -> TokenStream {
    let n_chi = res.n_chi;
    let n_ph = res.n_polar_h;
    let has_atoms = !res.atoms.is_empty();

    let pk = |s: &str| {
        if has_atoms {
            format_ident!("{s}")
        } else {
            format_ident!("_{s}")
        }
    };
    let pn = pk("n");
    let pca = pk("ca");
    let pc = pk("c");

    let chi_param = (n_chi > 0).then(|| quote! { , chi: [f32; #n_chi] });
    let ph_param = (n_ph > 0).then(|| quote! { , polar_h: [f32; #n_ph] });

    let imports = if !has_atoms {
        quote! {}
    } else if n_chi > 0 || n_ph > 0 {
        quote! {
            use crate::math::{TrigPair, sincosf};
            use crate::nerf::place;
        }
    } else {
        quote! {
            use crate::math::TrigPair;
            use crate::nerf::place;
        }
    };

    let chi_precompute: Vec<TokenStream> = (0..n_chi)
        .map(|i| {
            let var = format_ident!("chi{i}");
            quote! { let #var = sincosf(chi[#i]); }
        })
        .collect();

    let ph_precompute: Vec<TokenStream> = (0..n_ph)
        .map(|i| {
            let var = format_ident!("ph{i}");
            quote! { let #var = sincosf(polar_h[#i]); }
        })
        .collect();

    let placements: Vec<TokenStream> = res.atoms.iter().map(emit_placement).collect();

    let field_names: Vec<Ident> = res.atoms.iter().map(|a| field_ident(a.name)).collect();

    quote! {
        impl #ty {
            #[inline(always)]
            pub fn build(
                #pn: Vec3, #pca: Vec3, #pc: Vec3
                #chi_param
                #ph_param
            ) -> #coords {
                #imports
                #(#chi_precompute)*
                #(#ph_precompute)*
                #(#placements)*
                #coords { #(#field_names,)* }
            }
        }
    }
}

fn emit_placement(atom: &AtomSpec) -> TokenStream {
    let var = field_ident(atom.name);
    let ra = field_ident(atom.refs[0]);
    let rb = field_ident(atom.refs[1]);
    let rc = field_ident(atom.refs[2]);
    let d = f32_tok(atom.bond_length);
    let theta = trig_pair_tok(atom.bond_angle);

    match atom.torsion {
        TorsionSrc::Fixed(deg) => {
            let phi = trig_pair_tok(deg);
            quote! { let #var = place(#ra, #rb, #rc, #d, #theta, #phi); }
        }
        TorsionSrc::Chi(i) => {
            let phi = format_ident!("chi{i}");
            quote! { let #var = place(#ra, #rb, #rc, #d, #theta, #phi); }
        }
        TorsionSrc::PolarH(i, 0.0) => {
            let phi = format_ident!("ph{i}");
            quote! { let #var = place(#ra, #rb, #rc, #d, #theta, #phi); }
        }
        TorsionSrc::PolarH(i, offset) => {
            let phi_var = format_ident!("phi_{var}");
            let ph = format_ident!("ph{i}");
            let rad = offset * PI / 180.0;
            let co = f32_tok(rad.cos());
            let so = f32_tok(rad.sin());
            quote! {
                let #phi_var = TrigPair {
                    cos: #ph.cos * #co - #ph.sin * #so,
                    sin: #ph.sin * #co + #ph.cos * #so,
                };
                let #var = place(#ra, #rb, #rc, #d, #theta, #phi_var);
            }
        }
    }
}

fn validate_refs(res: &ResidueSpec) {
    let mut known: std::collections::HashSet<&str> = ["N", "CA", "C"].iter().copied().collect();

    for atom in res.atoms {
        for r in &atom.refs {
            assert!(
                known.contains(r),
                "{}: atom '{}' references '{}' which has not been placed yet. Known atoms: {:?}",
                res.name,
                atom.name,
                r,
                known
            );
        }

        match atom.torsion {
            TorsionSrc::Chi(i) => assert!(
                i < res.n_chi,
                "{}: atom '{}' uses Chi({}) but n_chi = {}",
                res.name,
                atom.name,
                i,
                res.n_chi
            ),
            TorsionSrc::PolarH(i, _) => assert!(
                i < res.n_polar_h,
                "{}: atom '{}' uses PolarH({}, _) but n_polar_h = {}",
                res.name,
                atom.name,
                i,
                res.n_polar_h
            ),
            TorsionSrc::Fixed(_) => {}
        }

        assert!(
            atom.bond_length > 0.8 && atom.bond_length < 2.5,
            "{}: atom '{}' has suspicious bond_length = {:.4} Å",
            res.name,
            atom.name,
            atom.bond_length
        );

        assert!(
            atom.bond_angle > 90.0 && atom.bond_angle < 180.0,
            "{}: atom '{}' has suspicious bond_angle = {:.4}°",
            res.name,
            atom.name,
            atom.bond_angle
        );

        known.insert(atom.name);
    }
}

fn emit_for_all_sidechains_macro() -> TokenStream {
    let arms: Vec<TokenStream> = ALL_RESIDUES
        .iter()
        .map(|res| {
            let ty = Ident::new(res.type_name, Span::call_site());
            let n_chi = Literal::usize_unsuffixed(res.n_chi);
            let n_ph = Literal::usize_unsuffixed(res.n_polar_h);
            let n = Literal::usize_unsuffixed(res.atoms.len());
            quote! { $callback!(#ty, #n_chi, #n_ph, #n); }
        })
        .collect();

    quote! {
        /// Invokes `$callback!(Type, N_CHI, N_POLAR_H, N)` for each of the 29
        /// sidechain types.
        ///
        /// The callback macro receives:
        /// - `Type` — the zero-sized residue type (e.g. `Ser`)
        /// - `N_CHI` — number of χ dihedral angles (literal)
        /// - `N_POLAR_H` — number of polar-hydrogen torsions (literal)
        /// - `N` — total number of sidechain atoms (literal)
        ///
        /// # Examples
        ///
        /// ```ignore
        /// macro_rules! print_name {
        ///     ($T:ident, $nc:literal, $np:literal, $n:literal) => {
        ///         println!("{}", <$T as Sidechain>::NAME);
        ///     };
        /// }
        /// for_all_sidechains!(print_name);
        /// ```
        #[macro_export]
        macro_rules! for_all_sidechains {
            ($callback:ident) => {
                #(#arms)*
            };
        }
    }
}
