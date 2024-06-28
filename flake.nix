{
  description = "Multi contact trajectory generation for the COM using Bezier curves";

  inputs = {
    nixpkgs.url = "github:nim65s/nixpkgs/gepetto";
    flake-parts = {
      url = "github:hercules-ci/flake-parts";
      inputs.nixpkgs-lib.follows = "nixpkgs";
    };
    hpp-centroidal-dynamics = {
      url = "github:humanoid-path-planner/hpp-centroidal-dynamics/release/5.1.0";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-parts.follows = "flake-parts";
    };
    ndcurves = {
      url = "github:loco-3d/ndcurves/release/1.5.0";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.flake-parts.follows = "flake-parts";
    };
  };

  outputs =
    inputs@{ flake-parts, ... }:
    flake-parts.lib.mkFlake { inherit inputs; } {
      imports = [ ];
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "aarch64-darwin"
        "x86_64-darwin"
      ];
      perSystem =
        {
          self',
          pkgs,
          system,
          ...
        }:
        {
          packages.default = pkgs.callPackage ./. {
            hpp-centroidal-dynamics = inputs.hpp-centroidal-dynamics.packages.${system}.default;
            ndcurves = inputs.ndcurves.packages.${system}.default;
          };
          devShells.default = pkgs.mkShell { inputsFrom = [ self'.packages.default ]; };
        };
    };
}
