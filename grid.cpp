#include "grid.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

Grid::Grid(int imax_init, int jmax_init, int boundary_size_init, double& PI, double& UI, double& VI, double& TI):
    _boundary_size(boundary_size_init),
    _imax_b(imax_init+2*boundary_size_init), // +2 for the total size of the vector / matrix
    _jmax_b(jmax_init+2*boundary_size_init),
    _imax(imax_init),
    _jmax(jmax_init),
    _fluid_cells_quantity(0)
{
    
    // Resizing the grid cells
    // the boundary size is given as a single value for all borders
    _cells.resize(Grid::imaxb(), std::vector<Cell>(Grid::jmaxb(), Cell(PI, UI, VI, TI)));

    // Resize Velocity Matrices
    _velocities[static_cast<int>(velocity_type::U)].resize(Grid::imaxb(), std::vector<double>(Grid::jmaxb(), UI));
    _velocities[static_cast<int>(velocity_type::V)].resize(Grid::imaxb(), std::vector<double>(Grid::jmaxb(), VI));

    // Setting the cells with a boundary cell
    // Top border
    std::vector<Cell> cells_temp;
    Grid::cells(cells_temp, matrix_selection::ROW, 0);
    std::for_each( cells_temp.begin(), cells_temp.end(), [](Cell& cell) { cell.set_border(border_position::BOTTOM); });
    Grid::set_cells(cells_temp, matrix_selection::ROW, 0);

    // Bottom
    Grid::cells(cells_temp, matrix_selection::ROW, jmaxb()-1);
    std::for_each( cells_temp.begin(), cells_temp.end(), [](Cell& cell) { cell.set_border(border_position::TOP); });
    Grid::set_cells(cells_temp, matrix_selection::ROW, jmaxb() - 1);

    // Left
    Grid::cells(cells_temp, matrix_selection::COLUMN, 0);
    std::for_each( cells_temp.begin(), cells_temp.end(), [](Cell& cell) { cell.set_border(border_position::LEFT); });
    Grid::set_cells(cells_temp, matrix_selection::COLUMN, 0);

    // RIGHT
    Grid::cells(cells_temp, matrix_selection::COLUMN, Grid::imaxb()-1);
    std::for_each( cells_temp.begin(), cells_temp.end(), [](Cell& cell) { cell.set_border(border_position::RIGHT); });
    Grid::set_cells(cells_temp, matrix_selection::COLUMN, Grid::imaxb() - 1);
};

int Grid::jmaxb() const {
    return _jmax_b;
};

int Grid::imaxb() const {
    return _imax_b;
};

int Grid::imax() const {
    return _imax;
};

int Grid::jmax() const {
    return _jmax;
};

Cell& Grid::cell(int i, int j) {
    return _cells.at(i).at(j);
};

void Grid::set_cell(Cell& cell, int i, int j) {
    _cells.at(i).at(j) = cell;
};

void Grid::velocity(matrix<double>& vec, velocity_type type) {
    if (!(type == velocity_type::U || type == velocity_type::V)) {
        std::cerr << "Wrong velocity type" << std::endl;
    }
    else {
        // Resize vector and set all values to 0
        vec.resize(Grid::imaxb(), std::vector<double>(Grid::jmaxb(), 0));
        for(int x=0; x < Grid::imaxb();x++) {
            for (int y=0; y < Grid::jmaxb(); y++) {
                // Accessing velocity
                vec.at(x).at(y) = _cells.at(x).at(y).velocity(type);
            }
        }
    }
}

void Grid::set_velocity(matrix<double>& vec, velocity_type type) {
    if (!(type == velocity_type::U || type == velocity_type::V)) {
        std::cerr << "Wrong velocity type" << std::endl;
    }
    else {
        for(int x=0; x< Grid::imaxb();x++) {
            for (int y=0; y < Grid::jmaxb(); y++) {
                // Setting velocity
                _cells.at(x).at(y).set_velocity(vec.at(x).at(y), type);
            }
        }
    }
}

void Grid::pressure(matrix<double>& vec) {

    // We don't need to do it every time, only once at first iteration. Some performance improvement.
    // static bool first_call_of_this_funtion = true; 

    // Resize vector and set all values to 0
    // if (first_call_of_this_funtion)
    vec.resize(Grid::imaxb(), std::vector<double>(Grid::jmaxb(), 0));

    // Iterate over cells
    for (int x = 0; x < Grid::imaxb(); x++) {
        for (int y = 0; y < Grid::jmaxb(); y++) {
            // Accessing pressure
            vec.at(x).at(y) = _cells.at(x).at(y).pressure();
        }
    }
    
    // first_call_of_this_funtion = false;
}

void Grid::set_pressure(matrix<double>& vec) {

    // Iterate over cells
    for(int x=0; x < Grid::imaxb();x++) {
        for (int y=0; y < Grid::jmaxb(); y++) {
            // Setting pressure
            _cells.at(x).at(y).set_pressure(vec.at(x).at(y));
        }
    }
}


void Grid::set_pressure_for_internal_boundaries() {

    static double pressure;
    pressure = 0; // Resetting this value to zero at every SOR iteration just to avoid issues

    // Iterate over internal cells
    for (int i = 1; i < Grid::imax() - 1; ++i) {
        for (int j = 1; j < Grid::jmax() - 1; ++j) {
            // Setting pressure for all boundary (NOSLIP) cells
            if (Grid::cell(i, j)._cellType == NOSLIP) {

                // Checking whether the neighbor cells belongs to FLUID, and taking its pressure

                // NORTH-EAST-ERN-cell
                if (((Grid::cell(i, j)._nbNorth)->_cellType == FLUID) &&
                    ((Grid::cell(i, j)._nbEast)->_cellType == FLUID))
                {
                    pressure = ((Grid::cell(i, j)._nbNorth)->pressure() + (Grid::cell(i, j)._nbEast)->pressure());
                    pressure /= 2;
                    Grid::cell(i, j).set_pressure(pressure);
                }

                // SOUTH-EAST-ERN-cell
                else if (((Grid::cell(i, j)._nbSouth)->_cellType == FLUID) &&
                    ((Grid::cell(i, j)._nbEast)->_cellType == FLUID))
                {
                    pressure = ((Grid::cell(i, j)._nbSouth)->pressure() + (Grid::cell(i, j)._nbEast)->pressure());
                    pressure /= 2;
                    Grid::cell(i, j).set_pressure(pressure);
                }

                // NORTH-WEST-ERN-cell
                else if (((Grid::cell(i, j)._nbNorth)->_cellType == FLUID) &&
                    ((Grid::cell(i, j)._nbWest)->_cellType == FLUID))
                {
                    pressure = ((Grid::cell(i, j)._nbNorth)->pressure() + (Grid::cell(i, j)._nbWest)->pressure());
                    pressure /= 2;
                    Grid::cell(i, j).set_pressure(pressure);
                }

                // SOUTH-WEST-ERN-cell
                else if (((Grid::cell(i, j)._nbSouth)->_cellType == FLUID) &&
                    ((Grid::cell(i, j)._nbWest)->_cellType == FLUID))
                {
                    pressure = ((Grid::cell(i, j)._nbSouth)->pressure() + (Grid::cell(i, j)._nbWest)->pressure());
                    pressure /= 2;
                    Grid::cell(i, j).set_pressure(pressure);
                }

                
                // NORTHERN-cell
                else if ((Grid::cell(i, j)._nbNorth)->_cellType == FLUID)
                    Grid::cell(i, j).set_pressure((Grid::cell(i, j)._nbNorth)->pressure());

                // EASTERN-cell
                else if ((Grid::cell(i, j)._nbEast)->_cellType == FLUID)
                    Grid::cell(i, j).set_pressure((Grid::cell(i, j)._nbEast)->pressure());

                // SOUTHERN-cell
                else if ((Grid::cell(i, j)._nbSouth)->_cellType == FLUID)
                    Grid::cell(i, j).set_pressure((Grid::cell(i, j)._nbSouth)->pressure());

                // WESTERN-cell (cowboy)
                else if ((Grid::cell(i, j)._nbWest)->_cellType == FLUID)
                    Grid::cell(i, j).set_pressure((Grid::cell(i, j)._nbWest)->pressure());
                

                /*
                // NORTHERN-cell
                else if ((Grid::cell(i, j)._nbNorth)->_cellType == FLUID)
                    Grid::cell(i, j).set_pressure(pressure);

                // EASTERN-cell
                else if ((Grid::cell(i, j)._nbEast)->_cellType == FLUID)
                    Grid::cell(i, j).set_pressure(pressure);

                // SOUTHERN-cell
                else if ((Grid::cell(i, j)._nbSouth)->_cellType == FLUID)
                    Grid::cell(i, j).set_pressure(pressure);

                // WESTERN-cell (cowboy)
                else if ((Grid::cell(i, j)._nbWest)->_cellType == FLUID)
                    Grid::cell(i, j).set_pressure(pressure);
                */

                // If there are no neighbor FLUID-cells, when taking the average pressure of four surrounding BOUNDARY cells
                else
                    Grid::cell(i, j).set_pressure(pressure);
            }
        }
    }
}

// Increment quantity of fluid cells by 1
void Grid::increment_fluid_cells() {
    _fluid_cells_quantity += 1;
};

// Get quantity of fluid cells
int Grid::get_fluid_cells_quantity() {
    return _fluid_cells_quantity;
}

void Grid::temperature(matrix<double>& vec) {
    // Resize vector and set all values to 0
    vec.resize(Grid::imaxb(), std::vector<double>(Grid::jmaxb(), 0));

    // Iterate over cells
    for (int x = 0; x < Grid::imaxb(); x++) {
        for (int y = 0; y < Grid::jmaxb(); y++) {
            // Accessing temperature
            vec.at(x).at(y) = _cells.at(x).at(y).temperature();
        }
    }
}

void Grid::set_temperature(matrix<double>& vec) {

    // Iterate over cells
    for (int x = 0; x < Grid::imaxb(); x++) {
        for (int y = 0; y < Grid::jmaxb(); y++) {
            // Setting temperature
            _cells.at(x).at(y).set_temperature(vec.at(x).at(y));
        }
    }
}

void Grid::cells(std::vector<Cell>& cells, matrix_selection m, int index) {
    // Resizing cells vector
    switch (m) {
    case matrix_selection::ROW:
        cells.resize(Grid::imaxb());
        for (int i = 0; i < Grid::imaxb(); i++) {
            cells.at(i) = Grid::cell(i, index);
        }
        break;
    case matrix_selection::COLUMN:
        cells.resize(Grid::jmaxb());
        for (int i = 0; i < Grid::jmaxb(); i++) {
            cells.at(i) = Grid::cell(index, i);
        }
        break;
    default:
        break;
    }
}

void Grid::innercells(std::vector<std::vector<Cell>> &cells) {
    cells.resize(Grid::imax(), std::vector<Cell>(Grid::jmax(), Cell()));
    for(int y = 0; y < Grid::jmax();y++) {
        for(int x = 0; x < Grid::imax();x++){
            // Accessing velocity
            cells.at(x).at(y) = Grid::_cells.at(x + Grid::_boundary_size).at(y + Grid::_boundary_size);
        }
    }
}

void Grid::set_innercells(std::vector<std::vector<Cell>> &cells) {
    for(int y = 0; y < Grid::jmax();y++) {
        for(int x = 0; x < Grid::imax();x++){
            // Accessing velocity
            Grid::_cells.at(x + Grid::_boundary_size).at(y + Grid::_boundary_size) = cells.at(x).at(y);
        }
    }
}

void Grid::set_cells(std::vector<Cell>& cells, matrix_selection m, int index) {
    // Resizing cells vector
    switch (m) {
    case matrix_selection::ROW:
        cells.resize(Grid::imaxb());
        for (int i = 0; i < Grid::imaxb(); i++) {
            Grid::set_cell(cells.at(i), i, index);
        }
        break;
    case matrix_selection::COLUMN:
        cells.resize(Grid::jmaxb());
        for (int i = 0; i < Grid::jmaxb(); i++) {
            Grid::set_cell(cells.at(i), index, i);
        }
        break;
    default:
        break;
    }
}

void Grid::print_velocity(velocity_type type) {
    if (!(type == velocity_type::U || type == velocity_type::V)) {
        std::cerr << "Wrong velocity type" << std::endl;
    }
    else {
        matrix<double> vel;
        Grid::velocity(vel, type);
        for(int y=Grid::jmaxb()-1;y >= 0;y--) {
            for(int x=0;x<Grid::imaxb();x++){
                // Accessing velocity
                std::cout<< vel.at(x).at(y) << " ";
            }
            // Print new line
            std::cout << std::endl;
        }
    }
}


void Grid::print_pressure() {
    for(int y=Grid::jmaxb()-1;y >= 0;y--) {
        for(int x=0;x<Grid::imaxb();x++){
                // Accessing pressure
                std::cout<< Grid::_cells.at(x).at(y).pressure() << " ";
            }
            // Print new line
            std::cout << std::endl;
        }
}

void Grid::print_temperature() {
    for (int y = Grid::jmaxb() - 1; y >= 0; y--) {
        for (int x = 0; x < Grid::imaxb(); x++) {
            // Accessing temperature
            std::cout << Grid::_cells.at(x).at(y).temperature() << " ";
        }
        // Print new line
        std::cout << std::endl;
    }
}


