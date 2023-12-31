// Solve the temperature equation

	fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian(alpha, T)
     ==
//        radiation->ST(rhoCpRef, T)
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

//    radiation->correct();

    fvOptions.correct(T);
