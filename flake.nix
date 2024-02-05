{
  description = "wtr-slider";

  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShell {
          packages = (with pkgs.rPackages; [
            pkgs.R
            languageserver
            httpgd
            pkgs.python3Packages.radian

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
          ++ (with pkgs; [ gettext ])
          ;
        };
      }
    );
}
