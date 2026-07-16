// Link libcint + libxc (built into ~/.local) and compile the libxc C shim.
use std::process::Command;

fn main() {
    let home = std::env::var("HOME").unwrap_or_else(|_| "/root".into());
    let inc = format!("{home}/.local/include");
    let lib = format!("{home}/.local/lib");
    let out = std::env::var("OUT_DIR").unwrap();
    let manifest = std::env::var("CARGO_MANIFEST_DIR").unwrap();

    // Compile the libxc shim (handles xc_func_type allocation in C).
    let shim_c = format!("{manifest}/src/xc_shim.c");
    let shim_o = format!("{out}/xc_shim.o");
    let status = Command::new("gcc")
        .args(["-c", "-O2", "-fPIC", &format!("-I{inc}"), &shim_c, "-o", &shim_o])
        .status()
        .expect("failed to run gcc for xc_shim.c");
    assert!(status.success(), "xc_shim.c compile failed");

    println!("cargo:rerun-if-changed={shim_c}");
    println!("cargo:rustc-link-search=native={lib}");
    println!("cargo:rustc-link-lib=dylib=cint");
    println!("cargo:rustc-link-lib=dylib=xc");
    println!("cargo:rustc-link-arg={shim_o}");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{lib}");
}
