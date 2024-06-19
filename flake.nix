{
  description = "wtr-slider";

  inputs.flake-parts.url = "github:hercules-ci/flake-parts";

  outputs = inputs@{ flake-parts, nixpkgs, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      perSystem = { self', pkgs, ... }:
        let
          overlay = final: prev: {
            rPackages = prev.rPackages.override {
              overrides = {
                unigd = prev.rPackages.unigd.overrideAttrs (attrs: {
                  env.NIX_LDFLAGS = "-lcairo";
                });
              };
            };
          };
          pkgs' = pkgs.extend overlay;
        in
        {
          devShells.default = pkgs'.mkShell {
            packages = (with pkgs'.rPackages; [
              pkgs'.R
              languageserver
              httpgd
              pkgs'.python3Packages.radian

              # dependencies for running the analysis code
              targets
              qs
              tidyverse
              jsonlite
              logr
              rstan
              bridgesampling
              brms
              bayestestR
            ])
            ++ (with pkgs'; [ gettext ]);
          };
        };
    };
}
