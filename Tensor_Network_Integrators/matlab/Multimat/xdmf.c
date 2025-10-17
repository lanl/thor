/*
 * XDMF file I/O
 */
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <hdf5.h>
#include <math.h>

#include "h5aux.h"
#include "xdmf.h"

void xdmf_new_file(const char *filename) {
    delete_file_if_exists(filename);
    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        perror("xdmf_new_file: Error creating file");
        return;
    }

    fprintf(file, "%s%s%s%s",
        "<Xdmf Version=\"3.0\">\n",
        "  <Domain><Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n",
        "  </Grid></Domain>\n",
        "</Xdmf>\n");

    fclose(file);
}


int count_lines_in_file(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("count_lines_in_file: error opening file");
        return -1;
    }

    int line_count = 0;
    char ch;

    while ((ch = fgetc(file)) != EOF) {
        if (ch == '\n') {
            line_count++;
        }
    }

    fclose(file);
    return line_count;
}


int read_file_into_lines(const char *filename, char ***lines) {
    // First, count the number of lines in the file
    int line_count = count_lines_in_file(filename);
    if (line_count < 0) {
        return -1;
    }

    // Allocate memory for the number of lines
    char **temp_lines = malloc(line_count * sizeof(char *));
    if (temp_lines == NULL) {
        perror("read_file_into_lines: error allocating memory");
        return -1;
    }

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("read_file_into_lines: error opening file");
        free(temp_lines);
        return -1;
    }

    char buffer[256];
    int current_line = 0;

    while (fgets(buffer, sizeof(buffer), file) != NULL) {
        // Strip newline character, if present
        size_t len = strlen(buffer);
        if (buffer[len - 1] == '\n') {
            buffer[len - 1] = '\0';
        }

        // Allocate memory for the line
        temp_lines[current_line] = malloc((len + 1) * sizeof(char));
        if (temp_lines[current_line] == NULL) {
            perror("read_file_into_lines: error allocating memory for line");
            fclose(file);
            for (int i = 0; i < current_line; i++) {
                free(temp_lines[i]);
            }
            free(temp_lines);
            return -1;
        }
        strcpy(temp_lines[current_line], buffer);
        current_line++;
    }

    fclose(file);

    // Assign the result to the output parameter
    *lines = temp_lines;
    return line_count;
}


FILE* xdmf_open_file_append(const char *filename) {
    FILE *file;
    int i, lines_count;
    char **lines;

    // read the file into lines
    lines_count = read_file_into_lines(filename, &lines);
    if (lines_count <= 0) {
        perror("xdmf_open_file_append: error reading lines");
        return NULL;  // Failure
    }

    // remove the file
    if (remove(filename)) {
        perror("delete_file_if_exists: error deleting file");
        return NULL;  // Failure
    }

    // create it again
    file = fopen(filename, "w");

    if (file == NULL) {
        perror("xdmf_open_file_append: error creating file");
        return NULL;
    }

    // write lines except the last two (closing)
    for (i=0; i < lines_count - 2; i++) {
        fprintf(file, "%s\n", lines[i]);
    }


    return file;
}

int xdmf_close_file(FILE* file) {
    if (file == NULL) {
        perror("xdmf_close_file: wrong file ID");
        return -1;
    }

    fprintf(file, "%s",
        "  </Grid></Domain>\n"
        "</Xdmf>\n");
    fclose(file);
    return 0;
}

void xdmf_write_mesh(FILE *file_id, const int dim, const double xl[],
        const double xr[], const int ncell[], const int nbdry) {
    double dx[3];
    int i;
    for (i=0; i<dim; ++i) {
        dx[i] = (xr[i] - xl[i])/(double)(ncell[i]);
    }
    if (dim == 2) {
        //  Number of cells: NX, NY; dimensions of the grid: (NY+1) (NX+1)"
        //  Origin: Y0, X0
        //  Spacing: DY, DX
        fprintf(file_id, "\n"
          "<Grid Name=\"StructuredGrid\">\n"
          " <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"%d %d\"/>\n"
          " <Geometry GeometryType=\"ORIGIN_DXDY\">\n"
          "  <DataItem DataType=\"Float\" Dimensions=\"2\" Format=\"XML\">\n"
          "    %14.7e %14.7e\n"
          "  </DataItem>\n"
          "  <DataItem DataType=\"Float\" Dimensions=\"2\" Format=\"XML\">\n"
          "    %14.7e %14.7e\n"
          "  </DataItem>\n"
          " </Geometry>\n",
          ncell[1] + 2*nbdry + 1, ncell[0] + 2*nbdry + 1,
          xl[1] - nbdry*dx[1], xl[0] - nbdry*dx[0], dx[1], dx[0]);
    }
    else if (dim == 3) {
        //  Number of cells: NX, NY, NZ;
        //  Dimensions of the grid: (NZ+1) (NY+1) (NX+1)"
        //  Origin: Z0, Y0, X0
        //  Spacing: DZ, DY, DX
        fprintf(file_id, "\n"
          "<Grid Name=\"StructuredGrid\">\n"
          " <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n"
          " <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
          "  <DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">\n"
          "    %14.7e %14.7e %14.7e\n"
          "  </DataItem>\n"
          "  <DataItem DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">\n"
          "    %14.7e %14.7e %14.7e\n"
          "  </DataItem>\n"
          " </Geometry>\n",
          ncell[2] + 2*nbdry + 1, ncell[1] + 2*nbdry + 1, ncell[0] + 2*nbdry + 1,
          xl[2] - nbdry*dx[2], xl[1] - nbdry*dx[1], xl[0] - nbdry*dx[0],
          dx[2], dx[1], dx[0]);
    }

}

void xdmf_write_dataset(FILE *file_id,
        const int dim, const int *dims,
        const double *xl, const double *xr,
        const char *fname_h5, const char *timestep_group,
        const char *varname, const int var_type) {
    double dx[3];
    int i;
    for (i=0; i<dim; ++i) {
        dx[i] = (xr[i] - xl[i])/(double)(dims[i]);
    }
    if (dim == 2) {
        if (var_type == CELL_VARIABLE) {
            //  Number of cells: NY, NX (+ 2*nbdry layers)
            fprintf(file_id, "\n"
              " <Attribute Name=\"%s\""
              "  AttributeType=\"Scalar\" Center=\"Cell\">\n"
              "   <DataItem Dimensions=\"%d %d\"\n"
              "     NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n"
              "     %s:/%s/%s\n"
              "   </DataItem>\n"
              " </Attribute>\n", varname,
              dims[1], dims[0],
              fname_h5, timestep_group, varname);
        }
        else if (var_type == NODE_VARIABLE) {
            fprintf(file_id, "\n"
              " <Attribute Name=\"%s\""
              "  AttributeType=\"Scalar\" Center=\"Node\">\n"
              "   <DataItem Dimensions=\"%d %d\"\n"
              "     NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n"
              "     %s:/%s/%s\n"
              "   </DataItem>\n"
              " </Attribute>\n", varname,
              dims[1], dims[0],
              fname_h5, timestep_group, varname);
        }
    }
    else if (dim == 3) {
        if (var_type == CELL_VARIABLE) {
            //  Number of cells: NZ, NY, NX (+ 2*nbdry layers)
            fprintf(file_id, "\n"
              " <Attribute Name=\"%s\""
              "  AttributeType=\"Scalar\" Center=\"Cell\">\n"
              "   <DataItem Dimensions=\"%d %d %d\"\n"
              "     NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n"
              "     %s:/%s/%s\n"
              "   </DataItem>\n"
              " </Attribute>\n", varname,
              dims[2], dims[1], dims[0],
              fname_h5, timestep_group, varname);
        }
        else if (var_type == NODE_VARIABLE) {
            fprintf(file_id, "\n"
              " <Attribute Name=\"%s\""
              "  AttributeType=\"Scalar\" Center=\"Node\">\n"
              "   <DataItem Dimensions=\"%d %d %d\"\n"
              "     NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n"
              "     %s:/%s/%s\n"
              "   </DataItem>\n"
              " </Attribute>\n", varname,
              dims[2], dims[1], dims[0],
              fname_h5, timestep_group, varname);
        }
    }
}

void xdmf_write_dataset_int(FILE *file_id,
        const int dim, const int *ncell,
        const double *xl, const double *xr,
        const char *fname_h5, const char *timestep_group,
        const char *varname, const int var_type) {
    double dx[3];
    int i;
    for (i=0; i<dim; ++i) {
        dx[i] = (xr[i] - xl[i])/(double)(ncell[i]);
    }
    if (dim == 2) {
        if (var_type == CELL_VARIABLE) {
            //  Number of cells: NY, NX (+ 2*nbdry layers)
            fprintf(file_id, "\n"
              " <Attribute Name=\"%s\""
              "  AttributeType=\"Scalar\" Center=\"Cell\">\n"
              "   <DataItem Dimensions=\"%d %d\"\n"
              "     NumberType=\"Integer\" Format=\"HDF\">\n"
              "     %s:/%s/%s\n"
              "   </DataItem>\n"
              " </Attribute>\n", varname,
              ncell[1], ncell[0],
              fname_h5, timestep_group, varname);
        }
        else if (var_type == NODE_VARIABLE) {
            fprintf(file_id, "\n"
              " <Attribute Name=\"%s\""
              "  AttributeType=\"Scalar\" Center=\"Node\">\n"
              "   <DataItem Dimensions=\"%d %d\"\n"
              "     NumberType=\"Integer\" Format=\"HDF\">\n"
              "     %s:/%s/%s\n"
              "   </DataItem>\n"
              " </Attribute>\n", varname,
              ncell[1] + 1, ncell[0] + 1,
              fname_h5, timestep_group, varname);
        }
    }
    else if (dim == 3) {
        if (var_type == CELL_VARIABLE) {
            //  Number of cells: NZ, NY, NX (+ 2*nbdry layers)
            fprintf(file_id, "\n"
              " <Attribute Name=\"%s\""
              "  AttributeType=\"Scalar\" Center=\"Cell\">\n"
              "   <DataItem Dimensions=\"%d %d %d\"\n"
              "     NumberType=\"Integer\" Format=\"HDF\">\n"
              "     %s:/%s/%s\n"
              "   </DataItem>\n"
              " </Attribute>\n", varname,
              ncell[2], ncell[1], ncell[0],
              fname_h5, timestep_group, varname);
        }
        else if (var_type == NODE_VARIABLE) {
            fprintf(file_id, "\n"
              " <Attribute Name=\"%s\""
              "  AttributeType=\"Scalar\" Center=\"Node\">\n"
              "   <DataItem Dimensions=\"%d %d %d\"\n"
              "     NumberType=\"Integer\" Format=\"HDF\">\n"
              "     %s:/%s/%s\n"
              "   </DataItem>\n"
              " </Attribute>\n", varname,
              ncell[2] + 1, ncell[1] + 1, ncell[0] + 1,
              fname_h5, timestep_group, varname);
        }
    }
}

void xdmf_write_timestamp(FILE *file_id, const double timestamp) {
    fprintf(file_id, "\n"
      " <Time Value=\"%14.7e\" />\n", timestamp);
}

void xdmf_close_group(FILE *file_id) {
    fprintf(file_id, "</Grid>\n\n");
}


