alphaRhoc = rho*alphac;
alphaRhocf = fvc::interpolate(alphaRhoc);
 alphaRhoPhic = alphacf*mixture.rhoPhi();

fvVectorMatrix UEqn
(

    fvm::ddt(alphaRhoc, U) + fvm::div(alphaRhoPhic, U)
  + turbulence->divDevTau(alphaRhoc,U)
  ==
     cloudSU
);
UEqn.relax();

fvOptions.constrain(UEqn);

volScalarField rAUc(1.0/UEqn.A());
volScalarField rASpUc(1.0/(UEqn.A() - cloudSUp));
surfaceScalarField rASpUcf("Dp", fvc::interpolate(rASpUc));

surfaceScalarField phicSUSu
(
    fvc::flux(rASpUc*cloudSUu)
);


surfaceScalarField phicSUSp
(
    fvc::interpolate(rASpUc*cloudSUp)
  );
if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
        (phicSUSu + phicSUSp*phi)/rASpUcf
      + (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
          - fvc::snGrad(p_rgh)
        )*mesh.magSf()
        )
      + (fvm::Sp(cloudSUp, U) - cloudSUp*U)
    );
    fvOptions.correct(U);
}
