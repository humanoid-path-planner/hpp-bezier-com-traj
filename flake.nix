{
  description = "Multi contact trajectory generation for the COM using Bezier curves";

  inputs.gepetto.url = "github:gepetto/nix";

  outputs =
    inputs:
    inputs.gepetto.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        overrideAttrs.hpp-bezier-com-traj = {
          src = lib.fileset.toSource {
            root = ./.;
            fileset = lib.fileset.unions [
              ./CMakeLists.txt
              ./include
              ./package.xml
              ./python
              ./src
              ./tests
            ];
          };
        };
      }
    );
}
