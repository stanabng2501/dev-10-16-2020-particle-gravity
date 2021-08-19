alphaRhoc = rho*alphac;
alphaRhocf = fvc::interpolate(alphaRhoc);
alphaRhoPhic = alphacf*mixture.rhoPhi();

fvVectorMatrix UcEqn
(

    fvm::ddt(alphaRhoc, Uc) + fvm::div(alphaRhoPhic, Uc)
  + turbulence->divDevTau(alphaRhoc,Uc)
  ==
     cloudSU
);

UcEqn.relax();

fvOptions.constrain(UcEqn);

volScalarField rAUc(1.0/UcEqn.A());
volScalarField rASpUc(1.0/(UcEqn.A() - cloudSUp));
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
        UcEqn
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
      + (fvm::Sp(cloudSUp, Uc) - cloudSUp*Uc)
    );
    fvOptions.correct(Uc);
}
