#include <iostream>
#include <cstring>
#include "core/mesh.hpp"

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [options] <input_file>" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --step1           Create nodes only" << std::endl;
    std::cout << "  --step2           Create nodes and simple triangulation" << std::endl;
    std::cout << "  --step3           Create nodes and Delaunay triangulation" << std::endl;
    std::cout << "  --step4           Create perturbed nodes and Delaunay triangulation" << std::endl;
    std::cout << "  --perturbation R  Set perturbation radius (default: 0.5)" << std::endl;
    std::cout << "  --help            Show this help message" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }
    
    int step = 3; // Default to step 3 (Delaunay triangulation)
    double perturbationRadius = 0.5;
    const char* inputFile = nullptr;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--help") == 0) {
            printUsage(argv[0]);
            return EXIT_SUCCESS;
        }
        else if (strcmp(argv[i], "--step1") == 0) {
            step = 1;
        }
        else if (strcmp(argv[i], "--step2") == 0) {
            step = 2;
        }
        else if (strcmp(argv[i], "--step3") == 0) {
            step = 3;
        }
        else if (strcmp(argv[i], "--step4") == 0) {
            step = 4;
        }
        else if (strcmp(argv[i], "--perturbation") == 0) {
            if (i + 1 < argc) {
                perturbationRadius = std::atof(argv[++i]);
            } else {
                std::cerr << "Error: --perturbation requires a value" << std::endl;
                return EXIT_FAILURE;
            }
        }
        else if (argv[i][0] != '-') {
            inputFile = argv[i];
        }
        else {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            printUsage(argv[0]);
            return EXIT_FAILURE;
        }
    }
    
    if (!inputFile) {
        std::cerr << "Error: Input file not specified" << std::endl;
        printUsage(argv[0]);
        return EXIT_FAILURE;
    }
    
    try {
        Mesh msh(inputFile);
        
        switch (step) {
            case 1:
                // Step 1: Create nodes only
                msh.createNodes();
                msh.printNodes();
                break;
                
            case 2:
                // Step 2: Create nodes and simple triangulation
                msh.createNodes();
                msh.printNodes();
                msh.simpleTriangulation();
                msh.printElements();
                break;
                
            case 3:
                // Step 3: Create nodes and Delaunay triangulation
                msh.createNodes();
                msh.delaunayTriangulation();
                msh.printNodes();
                msh.printElements();
                break;
                
            case 4:
                // Step 4: Create perturbed nodes and Delaunay triangulation
                msh.createNodes(perturbationRadius);
                msh.delaunayTriangulation();
                msh.printNodes();
                msh.printElements();
                break;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
