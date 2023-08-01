{
    // Scalar transport equation

    fvScalarMatrix sEqn
    (
        fvm::div(phi, s)
      - fvm::laplacian(Ds, s)
     ==
        fvOptions(s)
    );

    sEqn.relax();

    fvOptions.constrain(sEqn);

    sEqn.solve();

    fvOptions.correct(s);
}
