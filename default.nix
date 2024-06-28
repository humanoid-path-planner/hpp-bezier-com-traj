{
  lib,
  cddlib,
  clp,
  cmake,
  glpk,
  hpp-centroidal-dynamics,
  ndcurves,
  python3Packages,
  qpoases,
}:

python3Packages.buildPythonPackage {
  pname = "hpp-bezier-com-traj";
  version = "5.0.0";
  pyproject = false;

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

  strictDeps = true;

  nativeBuildInputs = [ cmake ];

  propagatedBuildInputs = [
    cddlib
    clp
    glpk
    hpp-centroidal-dynamics
    ndcurves
    qpoases
  ];

  cmakeFlags = [ "-DUSE_GLPK=ON" ];

  doCheck = true;

  pythonImportsCheck = [ "hpp_bezier_com_traj" ];

  meta = {
    description = "Multi contact trajectory generation for the COM using Bezier curves";
    homepage = "https://github.com/humanoid-path-planner/hpp-bezier-com-traj";
    license = lib.licenses.bsd2;
    maintainers = [ lib.maintainers.nim65s ];
  };
}
