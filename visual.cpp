#include "helper.hpp"
#include "visual.hpp"
#include <stdio.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTuple.h>
#include <vtkStructuredGridWriter.h>
#include "vector"

void write_vtkFile(
    const char *szProblem,
	int    timeStepNumber,
	double xlength,
    double ylength,
    int    imax,
    int    jmax,
	double dx,
	double dy,
    matrix <double> &U,
    matrix <double> &V,
    matrix <double> &P,
    matrix<double> &T)
{
    int i,j;
    char szFileName[80];
    FILE *fp=NULL;
    sprintf( szFileName, "%s.%i.vtk", szProblem, timeStepNumber );
    fp = fopen( szFileName, "w");
    if( fp == NULL )		       
    {
    char szBuff[80];
    sprintf( szBuff, "Failed to open %s", szFileName );
    ERROR( szBuff );
    return;
    }

    write_vtkHeader( fp, imax, jmax, dx, dy);
    write_vtkPointCoordinates(fp, imax, jmax, dx, dy);

    fprintf(fp,"POINT_DATA %i \n", (imax+1)*(jmax+1) );
	
    fprintf(fp,"\n");
    fprintf(fp, "VECTORS velocity float\n");
    for (j = 0; j < jmax + 1; j++)
    {
        for (i = 0; i < imax + 1; i++)
        {
            fprintf(fp, "%f %f 0\n", (U[i][j] + U[i][j + 1]) * 0.5, (V[i][j] + V[i + 1][j]) * 0.5);
        }
    }


  
    fprintf(fp,"\n");
    fprintf(fp,"CELL_DATA %i \n", ((imax)*(jmax)) );
    fprintf(fp, "SCALARS pressure float 1 \n"); 
    fprintf(fp, "LOOKUP_TABLE default \n");
    for(j = 1; j < jmax+1; j++)
    {
        for(i = 1; i < imax+1; i++)
        {
            fprintf(fp, "%f\n", P[i][j] );
        }
    }
  

  // TEMPERATURE
    fprintf(fp,"\n");
    fprintf(fp,"CELL_DATA %i \n", ((imax)*(jmax)) );
    fprintf(fp, "SCALARS temperature float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for(j = 1; j < jmax+1; j++)
    {
        for(i = 1; i < imax+1; i++)
        {
            fprintf(fp, "%f\n", T[i][j] );
        }
    }
    

    if( fclose(fp) )
    {
        char szBuff[80];
        sprintf( szBuff, "Failed to close %s", szFileName );
        ERROR( szBuff );
    }
}


void write_vtkHeader( FILE *fp, int imax, int jmax, 
                      double dx, double dy)
{
    if( fp == NULL )		       
    {
        char szBuff[80];
        sprintf( szBuff, "Null pointer in write_vtkHeader" );
        ERROR( szBuff );
        return;
    }

    fprintf(fp,"# vtk DataFile Version 2.0\n");
    fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");	
    fprintf(fp,"DATASET STRUCTURED_GRID\n");
    fprintf(fp,"DIMENSIONS  %i %i 1 \n", imax+1, jmax+1);
    fprintf(fp,"POINTS %i float\n", (imax+1)*(jmax+1) );
    fprintf(fp,"\n");
}


void write_vtkPointCoordinates(FILE *fp, int imax, int jmax, double dx, double dy)
{
    double originX = 0.0;  
    double originY = 0.0;

    int i = 0;
    int j = 0;

    for(j = 0; j < jmax+1; j++)
    {
        for(i = 0; i < imax+1; i++)
        {
            fprintf(fp, "%f %f 0\n", originX+(i*dx), originY+(j*dy) );
        }
    }
}

void VTKHelper::printVTKFile(
    Grid& grid,
    double dx,
    double dy,
    std::string casename,
    std::string outputdir,
    int timestep)
{

    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Lower Left Corner
    double x(0), y(0), z(0);
    for (int col = 0; col < grid.jmaxb() - 1; col++)
    {
        // Reset x again
        x = 0;
        for (int row = 0; row < grid.imaxb() - 1; row++)
        {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(grid.imaxb() - 1, grid.jmaxb() - 1, 1);
    structuredGrid->SetPoints(points);


    
    // Pressure Array
    vtkDoubleArray* Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);
    
    
    
    // TEMPERATURE
    // Temperature Array
    vtkDoubleArray* Temperature = vtkDoubleArray::New();
    Temperature->SetName("temperature");
    Temperature->SetNumberOfComponents(1);

    // DEPOSITED MATERIAL THICKNESS
    // Deposited Material Array
    vtkDoubleArray* DepositedMaterial = vtkDoubleArray::New();
    DepositedMaterial->SetName("deposition_thckns");
    DepositedMaterial->SetNumberOfComponents(1);

    // CONCENTRATION OF COMPONENT A
    // Concentration A Array
    vtkDoubleArray* ConcentrationA = vtkDoubleArray::New();
    ConcentrationA->SetName("concentrationA");
    ConcentrationA->SetNumberOfComponents(1);

    // CONCENTRATION OF COMPONENT B
    // Concentration B Array
    vtkDoubleArray* ConcentrationB = vtkDoubleArray::New();
    ConcentrationB->SetName("concentrationB");
    ConcentrationB->SetNumberOfComponents(1);

    // CONCENTRATION OF COMPONENT C
    // Concentration C Array
    vtkDoubleArray* ConcentrationC = vtkDoubleArray::New();
    ConcentrationC->SetName("concentrationC");
    ConcentrationC->SetNumberOfComponents(1);

    // CONCENTRATION OF COMPONENT D
    // Concentration D Array
    vtkDoubleArray* ConcentrationD = vtkDoubleArray::New();
    ConcentrationD->SetName("concentrationD");
    ConcentrationD->SetNumberOfComponents(1);


    // Velocity Array
    vtkDoubleArray* Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Geometry Array
    vtkDoubleArray* Geometry = vtkDoubleArray::New();
    Geometry->SetName("geometry");
    Geometry->SetNumberOfComponents(1);

    // Set number of tuples
    std::vector<std::vector<Cell>> cells;
    
    // PRESSURE
    std::vector<std::vector<double>> pressure;
    grid.pressure(pressure);

    // TEMPERATURE
    std::vector<std::vector<double>> temperature;
    grid.temperature(temperature);

    // DEPOSITION THICKNESS
    std::vector<std::vector<double>> deposition_thickness;
    grid.thickness(deposition_thickness);
    
    // CONCENTRATION of Component A
    std::vector<std::vector<double>> concentrationA;
    grid.concentration(concentrationA, ID::A);

    // CONCENTRATION of Component B
    std::vector<std::vector<double>> concentrationB;
    grid.concentration(concentrationB, ID::B);

    // CONCENTRATION of Component C
    std::vector<std::vector<double>> concentrationC;
    grid.concentration(concentrationC, ID::C);

    // CONCENTRATION of Component D
    std::vector<std::vector<double>> concentrationD;
    grid.concentration(concentrationD, ID::D);

    // TODO: create new matrix "geometry" depending on the type of cell and write it to the vtk-file
    // GEOMETRY OF OBSTACLES
    std::vector<std::vector<double>> geometry;
    geometry.resize(grid.imaxb(), std::vector<double>(grid.jmaxb(), 0));
    // write cell in formation to matrix
    for(int j = 0; j < grid.jmaxb(); j++) {
        for(int i = 0; i < grid.imaxb(); i++) {
            // print ones for obstacle cells
            if(grid.cell(i,j)._cellType == NOSLIP)
                geometry.at(i).at(j) = 1;
            // print ones for CATALYST cells
            else if (grid.cell(i, j)._cellType == CATALYST)
                geometry.at(i).at(j) = 2;
            // print zeros for non-obstacle cells
            else
                geometry.at(i).at(j) = 0;
        }
    }

    
    // Print pressure from bottom to top
    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            Pressure->InsertNextTuple(&pressure.at(i).at(j));
        }
    }
    

    // TEMPERATURE
    // Print temperature from bottom to top
    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            Temperature->InsertNextTuple(&temperature.at(i).at(j));
        }
    }

    // DEPOSITED MATERIAL THICKNESS
    // Print thickness from bottom to top
    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            DepositedMaterial->InsertNextTuple(&deposition_thickness.at(i).at(j));
        }
    }

    // CONCENTRATION of Component A
    // Print concentration from bottom to top
    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            ConcentrationA->InsertNextTuple(&concentrationA.at(i).at(j));
        }
    }

    // CONCENTRATION of Component B
    // Print concentration from bottom to top
    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            ConcentrationB->InsertNextTuple(&concentrationB.at(i).at(j));
        }
    }

    // CONCENTRATION of Component C
    // Print concentration from bottom to top
    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            ConcentrationC->InsertNextTuple(&concentrationC.at(i).at(j));
        }
    }

    // CONCENTRATION of Component D
    // Print concentration from bottom to top
    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            ConcentrationD->InsertNextTuple(&concentrationD.at(i).at(j));
        }
    }


    // GEOMETRY
    // Print geometry information from bottom to top
    for (int j = 1; j < grid.jmaxb() - 1; j++) {
        for (int i = 1; i < grid.imaxb() - 1; i++) {
            Geometry->InsertNextTuple(&geometry.at(i).at(j));
        }
    }





    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Get Velocity
    std::vector<std::vector<double>> velocity_U;
    std::vector<std::vector<double>> velocity_V;
    grid.velocity(velocity_U, velocity_type::U);
    grid.velocity(velocity_V, velocity_type::V);

    // Print Velocity from bottom to top
    for(int j = 1; j < grid.jmaxb(); j++)
    {
        for(int i = 1; i < grid.imaxb(); i++)
        {
            vel[0] = (velocity_U.at(i - 1).at(j) + velocity_U.at(i).at(j)) * 0.5;
            vel[1] = (velocity_V.at(i).at(j - 1) + velocity_V.at(i).at(j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    
    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);
    

    // TEMPERATURE
    // Add Temperature to Structured Grid
    structuredGrid->GetCellData()->AddArray(Temperature);

    // DEPOSITED MATERIAL THICKNESS
    // Add Thickness to Structured Grid
    structuredGrid->GetCellData()->AddArray(DepositedMaterial);
    
    // CONCENTRATION OF COMPONENT A
    // Add Concentration to Structured Grid
    structuredGrid->GetCellData()->AddArray(ConcentrationA);

    // CONCENTRATION OF COMPONENT B
    // Add Concentration to Structured Grid
    structuredGrid->GetCellData()->AddArray(ConcentrationB);

    // CONCENTRATION OF COMPONENT C
    // Add Concentration to Structured Grid
    structuredGrid->GetCellData()->AddArray(ConcentrationC);

    // CONCENTRATION OF COMPONENT D
    // Add Concentration to Structured Grid
    structuredGrid->GetCellData()->AddArray(ConcentrationD);

    // Add Geometry to Structured Grid
    structuredGrid->GetCellData()->AddArray(Geometry);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);



    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname = outputdir + '/' + casename + "." + std::to_string(timestep) + ".vtk";
    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}