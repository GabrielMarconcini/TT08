{
    // Second specie mass fraction transport

    fvScalarMatrix Y2Eqn
    (
        fvm::div(phi, Y2)
      - fvm::laplacian(D, Y2)
     ==
        fvOptions(Y2)
    );

    Y2Eqn.relax();

    fvOptions.constrain(Y2Eqn);

    Y2Eqn.solve();

    fvOptions.correct(Y2);
}
