function envActivate(activOption, envPath, requiredPackages)
    # Function that creates and activates an environment in which to run the Package
    # activeOption=1 individually checks which required packages are already installed
    # activOption=2 installs all required packages regardless of which ones already exist in the environment

    # Activate the environment
    Pkg.activate(envPath)
    
    if activOption == 1

        ## PACKAGE INSTALLATION OPTION 1

        # Get the current environment's package names
        installed_packages = keys(Pkg.installed());
        # Flag to track if any package is missing
        missing_packages = false;
        # Check for each required package and add if not already in Project.toml
        for pkg in requiredPackages
            if !(pkg in installed_packages)
                global missing_packages = true;
                println("Missing package: $pkg. Installing...")  # Print the missing package name
                Pkg.add(pkg);
            elseif pkg == "StructuralIdentifiability" 
                println("Specific version of $pkg required. Installing...")  # Print the package name
                Pkg.add(Pkg.PackageSpec(;name="StructuralIdentifiability", version="v0.5.1")); # need to use older version
            elseif pkg == "CSV"
                println("Specific version of $pkg required. Installing...")  # Print the package name
                Pkg.add(pkg);
            elseif pkg == "Random"
                println("Specific version of $pkg required. Installing...")  # Print the package name
                Pkg.add(Pkg.PackageSpec(;name="Random", version="1.9.3"));
            else
                println("$pkg is already installed.")  # Print the already installed package
            end
        end
    elseif activOption == 2

        ## PACKAGE INSTALLATION OPTION 2

        # Ensure that the required packages are installed
        Pkg.add(requiredPackages);
        Pkg.add(Pkg.PackageSpec(;name="StructuralIdentifiability", version="v0.5.1")); # need to use older version
        Pkg.add(Pkg.PackageSpec(;name="Random", version="1.9.3")); # need to use older version
    end

    # Check that Manifest.toml is consistent with Project.tomls of the current project and all its
    # dependencies, updating it if necessary, then instantiate. 
    Pkg.resolve();

end
    