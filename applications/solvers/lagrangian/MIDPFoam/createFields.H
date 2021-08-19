#include "readGravitationalAcceleration.H"


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
#include "readhRef.H"
#include "gh.H"


Info<< "Reading/calculating continuous-phase face flux field phi\n"
    << endl;
    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    );

Info << "Creating field alphac\n" << endl;

volScalarField alphac
(
    IOobject
    (
        "alphac",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 0)
);

Info<< "Creating turbulence model\n" << endl;


multiphaseEmulsionMix mixture(U, phi);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alphac*mixture.rho()
);
rho.oldTime();





Info<< "Constructing clouds" << endl;
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alphac*mixture.mu()
);

parcelCloudList clouds(rho, U, mu, g);

scalar alphacMin
(
    1 - mesh.solverDict(alphac.name()).lookup<scalar>("max")
);

alphac = max(1.0 - clouds.theta(), alphacMin);
alphac.correctBoundaryConditions();

surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));

surfaceScalarField alphaPhic
(
    "alphaPhic",
    alphacf*phi
);


volScalarField alphaRhoc
(    IOobject
    (
        "alphaRhoc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*alphac
);


surfaceScalarField alphaRhocf
(
    IOobject
    (
        "alphaRhocf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(alphaRhoc)
);

surfaceScalarField alphaRhoPhic
(
  "alphaRhoPhic",
  alphacf*mixture.rhoPhi()
);
autoPtr<incompressible::momentumTransportModel> turbulence
(
incompressible::momentumTransportModel::New(U, phi, mixture)
);




surfaceScalarField rhocf
(
    IOobject
    (
        "rhocf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(rho)
);


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);



label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
   p_rgh = p - rho*gh;

}


mesh.setFluxRequired(p_rgh.name());


#include "createFvOptions.H"
