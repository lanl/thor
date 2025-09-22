/*
 * HDF5 file I/O
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <hdf5.h>
#include <math.h>

// Error-checking macro
#define CHECK_H5ERROR(call) do { \
        herr_t status = (call);\
        if (status < 0) {      \
            fprintf(stderr,    \
            "Error at line %d in function %s\n", \
             __LINE__, __func__); \
            H5Eprint(H5E_DEFAULT, stderr); \
            exit(EXIT_FAILURE);\
        }} while(0)

int delete_file_if_exists(const char *filename) {
    if (access(filename, F_OK) == 0) {
        if (remove(filename) == 0)
            return 0;  // Success
        else {
            perror("delete_file_if_exists: error deleting file");
            return -1;  // Failure
        }
    }
    return 1;  // File not found
}


hid_t h5_new_file(const char *filename) {
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_H5ERROR(file_id);
    return file_id;
}

hid_t h5_open_existing_rdwr(const char *filename) {
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    CHECK_H5ERROR(file_id);
    return file_id;
}


hid_t h5_open_existing_rdonly(const char *filename) {
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    CHECK_H5ERROR(file_id);
    return file_id;
}


herr_t count_groups(hid_t group_id, const char *name, const H5L_info_t *info, void *op_data) {
    int *group_count = (int *)op_data;

    // Check if the object is a group
    H5O_info_t object_info;
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR >= 14
    // Version 1.14 or higher (5 arguments for H5Oget_info_by_name)
    CHECK_H5ERROR(H5Oget_info_by_name(group_id, name, &object_info, H5O_INFO_BASIC, H5P_DEFAULT));
#else
    // Version lower than 1.14 (4 arguments for H5Oget_info_by_name)
    CHECK_H5ERROR(H5Oget_info_by_name(group_id, name, &object_info, H5P_DEFAULT));
#endif    

    if (object_info.type == H5O_TYPE_GROUP) {
        (*group_count)++;
    }

    return 0;  // Continue iteration
}


int h5_count_groups(const hid_t fogr_id) {

    int num_groups;

    num_groups = 0;
    CHECK_H5ERROR(H5Literate(fogr_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, count_groups, &num_groups));

    return num_groups;

}


hid_t h5_open_group(hid_t file_id, const char *group_name) {
    hid_t group_id;

    // Check if the group already exists
    if (H5Lexists(file_id, group_name, H5P_DEFAULT) > 0) {
        // If the group exists, open it
        group_id = H5Gopen(file_id, group_name, H5P_DEFAULT);
        CHECK_H5ERROR(group_id);
    } else {
        // If the group does not exist, create it
        group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        CHECK_H5ERROR(group_id);
    }

    // Return the group identifier
    return group_id;
}

// Function to write a 1D array to an HDF5 file
herr_t h5_write_1d(hid_t fogr_id, const char *dataset_name,
                   double *data, const int len) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[1] = {(hsize_t)len};

    // Define the dimensions of the 1D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_1d: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_DOUBLE, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_1d: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    if (status < 0) {
        perror("h5_write_1d: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}


// Function to write a 1D array of integers to an HDF5 file
herr_t h5_write_1d_int(hid_t fogr_id, const char *dataset_name,
                   int *data, const int len) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[1] = {(hsize_t)len};

    // Define the dimensions of the 2D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_1d_int: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_INT, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_1d_int: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0]);
    if (status < 0) {
        perror("h5_write_1d_int: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}

// Function to write a 2D array to an HDF5 file
herr_t h5_write_2d(hid_t fogr_id, const char *dataset_name,
                   double **data, const int ncell[]) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[2] = {(hsize_t)ncell[0], (hsize_t)ncell[1]};

    // Define the dimensions of the 2D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_2d: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_DOUBLE, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_2d: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0][0]);
    if (status < 0) {
        perror("h5_write_2d: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}


// Function to write a 2D array of integers to an HDF5 file
herr_t h5_write_2d_int(hid_t fogr_id, const char *dataset_name,
                   int **data, const int ncell[]) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[2] = {(hsize_t)ncell[0], (hsize_t)ncell[1]};

    // Define the dimensions of the 2D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_2d_int: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_INT, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_2d_int: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0][0]);
    if (status < 0) {
        perror("h5_write_2d_int: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}


// Function to write a 2D array of long long to an HDF5 file
herr_t h5_write_2d_llong(hid_t fogr_id, const char *dataset_name,
                   long long **data, const int ncell[]) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[2] = {(hsize_t)ncell[0], (hsize_t)ncell[1]};

    // Define the dimensions of the 2D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_2d_llong: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_LLONG, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_2d_llong: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0][0]);
    if (status < 0) {
        perror("h5_write_2d_int: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}


// Function to write a 3D array to an HDF5 file
herr_t h5_write_3d(hid_t fogr_id, const char *dataset_name,
                   double ***data, const int ncell[]) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[3] = {(hsize_t)ncell[0],
                       (hsize_t)ncell[1],
                       (hsize_t)ncell[2]};

    // Define the dimensions of the 2D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_3d: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_DOUBLE, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_3d: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0][0][0]);
    if (status < 0) {
        perror("h5_write_3d: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}


// Function to write a 3D integer array to an HDF5 file
herr_t h5_write_3d_int(hid_t fogr_id, const char *dataset_name,
                   int ***data, const int ncell[]) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[3] = {(hsize_t)ncell[0],
                       (hsize_t)ncell[1],
                       (hsize_t)ncell[2]};

    // Define the dimensions of the 2D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_3d_int: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_INT, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_3d_int: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0][0][0]);
    if (status < 0) {
        perror("h5_write_3d_int: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}


// Function to write a 3D long long integer array to an HDF5 file
herr_t h5_write_3d_llong(hid_t fogr_id, const char *dataset_name,
                   long long ***data, const int ncell[]) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[3] = {(hsize_t)ncell[0],
                       (hsize_t)ncell[1],
                       (hsize_t)ncell[2]};

    // Define the dimensions of the 2D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_3d_llong: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_LLONG, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_3d_llong: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0][0][0]);
    if (status < 0) {
        perror("h5_write_3d_llong: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}


// Function to write a 4D array to an HDF5 file
herr_t h5_write_4d(hid_t fogr_id, const char *dataset_name,
                   double ****data, const int ncell[]) {

    // Convert from int to hsize_t
    int i;
    hsize_t dims[4] = {(hsize_t)ncell[0],
                       (hsize_t)ncell[1],
                       (hsize_t)ncell[2],
                       (hsize_t)ncell[3]};

    // Define the dimensions of the 4D array (dataspace)
    hid_t dataspace_id = H5Screate_simple(4, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_4d: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataspace
    }

    // Create the dataset
    hid_t dataset_id = H5Dcreate(fogr_id, dataset_name,
                                 H5T_NATIVE_DOUBLE, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        perror("h5_write_4d: Error creating dataset");
        H5Sclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    // Write the 2D array data to the dataset
    herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE,
                             H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0][0][0][0]);
    if (status < 0) {
        perror("h5_write_4d: Error writing data to dataset");
    }

    // Close the dataset, dataspace, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    return status;  // Return success (0) or failure (-1)
}


double** h5_read_2d(hid_t group_id, const char *dataset_name, int *dims_exp) {
    hid_t dataset_id, dataspace_id, memspace_id;
    hsize_t dims[2];  // For storing the dimensions of the dataset
    double **array;

    // Open the dataset
    dataset_id = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error in h5_read_2d: Could not open dataset %s\n", dataset_name);
        exit(EXIT_FAILURE);
    }

    // Get the dataset's dataspace and check its dimensions
    dataspace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

    // Check if dimensions match the expected size (dims_exp[0] x dims_exp[1])
    if (dims[0] != dims_exp[0] || dims[1] != dims_exp[1]) {
        fprintf(stderr, "Error in h5_read_2d: Dataset %s has incorrect dimensions"
            " (expected %d x %d, got %llu x %llu)\n",
            dataset_name, dims_exp[0], dims_exp[1], dims[0], dims[1]);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        exit(EXIT_FAILURE);
    }

    // Allocate array
    array = (double **) malloc(dims[1] * sizeof(double *));
    for (int j = 0; j < dims[1]; j++) {
        if (j == 0)
            array[0] = (double *) malloc(dims[0]*dims[1] * sizeof(double));
        else
            array[j]  =  array[j-1] + dims[0];
    }

    // Create a memory space for the dataset
    hsize_t mem_dims[2] = {dims[0], dims[1]};
    memspace_id = H5Screate_simple(2, mem_dims, NULL);

    // Read the dataset into the preallocated 2D array
    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, &(array[0][0]));

    // Close resources
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    return array;
} // h5_read_2d



int** h5_read_2d_int(hid_t group_id, const char *dataset_name, int *dims_exp) {
    hid_t dataset_id, dataspace_id, memspace_id;
    hsize_t dims[2];  // For storing the dimensions of the dataset
    int **array;

    // Open the dataset
    dataset_id = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
    if (dataset_id < 0) {
        fprintf(stderr, "Error in h5_read_2d_int: Could not open dataset %s\n", dataset_name);
        exit(EXIT_FAILURE);
    }

    // Get the dataset's dataspace and check its dimensions
    dataspace_id = H5Dget_space(dataset_id);
    H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

    // Check if dimensions match the expected size (dims_exp[0] x dims_exp[1])
    if (dims[0] != dims_exp[0] || dims[1] != dims_exp[1]) {
        fprintf(stderr, "Error in h5_read_2d_int: Dataset %s has incorrect dimensions"
            " (expected %d x %d, got %llu x %llu)\n",
            dataset_name, dims_exp[0], dims_exp[1], dims[0], dims[1]);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        exit(EXIT_FAILURE);
    }

    // Allocate array
    array = (int **) malloc(dims[1] * sizeof(int *));
    for (int j = 0; j < dims[1]; j++) {
        if (j == 0)
            array[0] = (int *) malloc(dims[0]*dims[1] * sizeof(int));
        else
            array[j]  =  array[j-1] + dims[0];
    }

    // Create a memory space for the dataset
    hsize_t mem_dims[2] = {dims[0], dims[1]};
    memspace_id = H5Screate_simple(2, mem_dims, NULL);

    // Read the dataset into the preallocated 2D array
    H5Dread(dataset_id, H5T_NATIVE_INT, memspace_id, dataspace_id, H5P_DEFAULT, &(array[0][0]));

    // Close resources
    H5Sclose(memspace_id);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);

    return array;
} // h5_read_2d_int



herr_t h5_write_attr(hid_t fogr_id, const char *attr_name, double val) {


    hid_t attr_id = H5Acreate(fogr_id, attr_name, H5T_NATIVE_DOUBLE, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id < 0) {
        perror("h5_write_attr: Error creating attribute");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &val);
    H5Aclose(attr_id);
    return 0;

}

herr_t h5_write_attr_int(hid_t fogr_id, const char *attr_name, int ival) {


    hid_t attr_id = H5Acreate(fogr_id, attr_name, H5T_NATIVE_INT, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id < 0) {
        perror("h5_write_attr_int: Error creating attribute");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }
    H5Awrite(attr_id, H5T_NATIVE_INT, &ival);
    H5Aclose(attr_id);
    return 0;

}

herr_t h5_write_attr_1d(hid_t fogr_id, const char *attr_name, double *data, const int len) {


    hsize_t dims[1] = {(hsize_t)len};
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_attr: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    hid_t attr_id = H5Acreate(fogr_id, attr_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id < 0) {
        perror("h5_write_attr: Error creating attribute");
        H5Fclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }
    H5Awrite(attr_id, H5T_NATIVE_DOUBLE, data);
    H5Aclose(attr_id);
    H5Sclose(dataspace_id);
    return 0;

}

herr_t h5_write_attr_1d_int(hid_t fogr_id, const char *attr_name, int *data, const int len) {


    hsize_t dims[1] = {(hsize_t)len};
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    if (dataspace_id < 0) {
        perror("h5_write_attr: Error creating dataspace");
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }

    hid_t attr_id = H5Acreate(fogr_id, attr_name, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attr_id < 0) {
        perror("h5_write_attr: Error creating attribute");
        H5Fclose(dataspace_id);
        H5Fclose(fogr_id);
        return -1;  // Error creating dataset
    }
    H5Awrite(attr_id, H5T_NATIVE_INT, data);
    H5Aclose(attr_id);
    H5Sclose(dataspace_id);
    return 0;

}


int h5_read_attr_int(hid_t group_id, const char *attr_name) {
    int attr_value;
    hid_t attr_id, space_id, attr_type;

    // Open the attribute
    attr_id = H5Aopen(group_id, attr_name, H5P_DEFAULT);
    if (attr_id < 0) {
        fprintf(stderr, "Error in h5_read_attr_int: could not open attribute %s\n", attr_name);
        exit(-1);
    }

    // Get attribute's datatype and dataspace
    attr_type = H5Aget_type(attr_id);
    space_id = H5Aget_space(attr_id);

    // Check if attribute is a scalar
    if (H5Sget_simple_extent_ndims(space_id) != 0) {
        fprintf(stderr, "Error: Attribute %s is not scalar\n", attr_name);
        exit(EXIT_FAILURE);
    }

    // Read the attribute
    H5Aread(attr_id, attr_type, &attr_value);

    // Close resources
    H5Sclose(space_id);
    H5Tclose(attr_type);
    H5Aclose(attr_id);

    return attr_value;
}



void h5_read_attr_1d(hid_t group_id, const char *attr_name, double *attr_data, const int dim) {
    hid_t attr_id, space_id, attr_type;
    hsize_t attr_dims[1];  // For storing the dimensions of the attribute

    // Open the attribute
    attr_id = H5Aopen(group_id, attr_name, H5P_DEFAULT);
    if (attr_id < 0) {
        fprintf(stderr, "Error in h5_read_attr_1d:"
            " Could not open attribute %s\n", attr_name);
        exit(EXIT_FAILURE);
    }

    // Get the attribute's dataspace and check its dimensions
    space_id = H5Aget_space(attr_id);
    H5Sget_simple_extent_dims(space_id, attr_dims, NULL);

    if (attr_dims[0] != dim) {
        fprintf(stderr, "Error in h5_read_attr_1d:"
            " Attribute %s has incorrect dimensions (expected %d, got %llu)\n",
            attr_name, dim, attr_dims[0]);
        H5Sclose(space_id);
        H5Aclose(attr_id);
        exit(EXIT_FAILURE);
    }

    // Get attribute's datatype
    attr_type = H5Aget_type(attr_id);

    // Read the attribute data into the preallocated array
    H5Aread(attr_id, attr_type, attr_data);

    // Close resources
    H5Tclose(attr_type);
    H5Sclose(space_id);
    H5Aclose(attr_id);
}


void h5_read_attr_1d_int(hid_t group_id, const char *attr_name, int *attr_data, const int dim) {
    hid_t attr_id, space_id, attr_type;
    hsize_t attr_dims[1];  // For storing the dimensions of the attribute

    // Open the attribute
    attr_id = H5Aopen(group_id, attr_name, H5P_DEFAULT);
    if (attr_id < 0) {
        fprintf(stderr, "Error in h5_read_attr_1d_int:"
            " Could not open attribute %s\n", attr_name);
        exit(EXIT_FAILURE);
    }

    // Get the attribute's dataspace and check its dimensions
    space_id = H5Aget_space(attr_id);
    H5Sget_simple_extent_dims(space_id, attr_dims, NULL);

    if (attr_dims[0] != dim) {
        fprintf(stderr, "Error in h5_read_attr_1d_int:"
            " Attribute %s has incorrect dimensions (expected %d, got %llu)\n",
            attr_name, dim, attr_dims[0]);
        H5Sclose(space_id);
        H5Aclose(attr_id);
        exit(EXIT_FAILURE);
    }

    // Get attribute's datatype
    attr_type = H5Aget_type(attr_id);

    // Read the attribute data into the preallocated array
    H5Aread(attr_id, attr_type, attr_data);

    // Close resources
    H5Tclose(attr_type);
    H5Sclose(space_id);
    H5Aclose(attr_id);
}


void copy_velocity_component_2d (double ***vel_for_2dnode,
        double ***var_for_2dnode, const int c, const int *dims) {
    int i, j;
    for (j=0; j<dims[1]; ++j) {
    for (i=0; i<dims[0]; ++i) {
            (*var_for_2dnode)[j][i] = vel_for_2dnode[j][i][c];
    }}
}


void copy_velocity_component_3d (double ****vel_for_3dnode,
        double ****var_for_3dnode, const int c, const int *dims) {
    int i, j, k;
    for (k=0; k<dims[2]; ++k) {
    for (j=0; j<dims[1]; ++j) {
    for (i=0; i<dims[0]; ++i) {
        (*var_for_3dnode)[k][j][i] = vel_for_3dnode[k][j][i][c];
    }}}
}



