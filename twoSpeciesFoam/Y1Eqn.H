{
    // First specie mass fraction transport

    fvScalarMatrix Y1Eqn
    (
        fvm::div(phi, Y1)
      - fvm::laplacian(D, Y1)
     ==
        fvOptions(Y1)
    );

    Y1Eqn.relax();

    fvOptions.constrain(Y1Eqn);

    Y1Eqn.solve();

    fvOptions.correct(Y1);
}
