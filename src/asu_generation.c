#include "asu_utils.h"
#include "asu.h"
#include "asu_generation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "mpi.h"
#include "algebra.h"
#include "molecule.h"
#include "molecule_utils.h"
#include "randomgen.h"
#include "check_structure.h"   // convert_atom2atom_vdw, find_mol_len

int asu_try_one_random_place(asu *unit, float box[3][3], molecule *mol)
{
       float random_rot[3][3];
       float random_trans[3];
       float inv_lat_vec[3][3];
       float lat_vec_trans[3][3];

       inverse_mat3b3(inv_lat_vec, box);
       mat3b3_transpose(inv_lat_vec, inv_lat_vec);
       mat3b3_transpose(lat_vec_trans, box);

       int at = 0;
       for(int m = 0; m < unit->n_mol_types; m++)
       {
              for(int st = 0; st < unit->stoic[m]; st++)
              {
                     generate_random_rotation_matrix(random_rot);
                     generate_random_translation_vector(random_trans);

                     for(int mat = 0; mat < unit->n_atoms_in_mol[m]; mat++)
                     {
                            float atom[3] = {mol[m].X[mat], mol[m].Y[mat], mol[m].Z[mat]};
                            vector3_mat3b3_multiply(random_rot, atom, atom);

                            // Convert to fractional and translate
                            vector3_mat3b3_multiply(inv_lat_vec, atom, atom);
                            vector3_add(atom, random_trans, atom);

                            vector3_mat3b3_multiply(lat_vec_trans, atom, atom);
                            unit->Xcord[at] = atom[0];
                            unit->Ycord[at] = atom[1];
                            unit->Zcord[at] = atom[2];
                            unit->atoms[2*at] = mol[m].atoms[2*mat];
                            unit->atoms[2*at + 1] = mol[m].atoms[2*mat + 1];
                            at++;
                     }
              }
       }
       return 1;
}

float cal_asu_two_mol_sr(asu *unit, int ith_mol, int jth_mol, float *atoms_vdw)
{
      int index_i = unit->mol_index[ith_mol];
      int index_j = unit->mol_index[jth_mol];
      float sr = 1000000.0;
      float temp = 0.0;
      float dx, dy, dz;
      for(int i = 0; i < unit->n_atoms_in_mol[unit->mol_types[ith_mol]]; i++)
      {
              for(int j = 0; j < unit->n_atoms_in_mol[unit->mol_types[jth_mol]]; j++)
              {
                     dx = unit->Xcord[index_i + i] - unit->Xcord[index_j + j];
                     dy = unit->Ycord[index_i + i] - unit->Ycord[index_j + j];
                     dz = unit->Zcord[index_i + i] - unit->Zcord[index_j + j];
                     temp = sqrt(dx * dx + dy * dy + dz * dz) / (atoms_vdw[index_i + i] + (atoms_vdw[index_j + j]));
                     if (temp < sr)
                     {
                         sr = temp;
                     }
              }

      }
      return sr;
}

// Smallest specific radius (min interatomic distance / vdW-radius sum) over
// all molecule pairs. Used to reject overlapping placements. Requires >= 2
// molecules; with a single molecule there is no pair and sr stays "infinite".
float cal_asu_sr(asu *unit)
{
       float atoms_vdw[unit->n_atoms];
       float sr = 1000000.0;
       float temp = 0;
       convert_atom2atom_vdw(unit->atoms, atoms_vdw, unit->n_atoms);

       for(int i = 0; i < unit->n_mols - 1; i++){

              for(int j = i + 1; j < unit->n_mols; j++)
              {
                     temp = cal_asu_two_mol_sr(unit, i, j, atoms_vdw);
                     if (temp < sr)
                     {
                            sr = temp;
                     }

              }
       }
       return sr;

}

void estimate_box(float box[3][3], int *stoic, molecule *mol, int mol_types){

     float mol_length;
     float estimate_vol = 0;
     float a;
     for(int i = 0; i < mol_types; i++)
     {
         mol_length = find_mol_len(mol[i].X, mol[i].Y, mol[i].Z, mol[i].num_of_atoms);
         estimate_vol += mol_length * mol_length * mol_length * stoic[i];

     }

     estimate_vol = estimate_vol * 1.5 * 2;
     a = cbrt(estimate_vol);
     box[0][0] = a;
     box[1][1] = a;
     box[2][2] = a;

     box[0][1] = 0;
     box[0][2] = 0;
     box[1][2] = 0;

     box[1][0] = 0;
     box[2][0] = 0;
     box[2][1] = 0;
}

int generate_one_asu(asu *unit, molecule *mol, float sr_min, float sr_max, int max_attempts)
// *unit is an initialized asymmetric unit
{
       float box[3][3];
       float sr;

       estimate_box(box, unit->stoic, mol, unit->n_mol_types);

       for(int i = 0; i < max_attempts; i++){
              asu_try_one_random_place(unit, box, mol);
              sr = cal_asu_sr(unit);

              if ((sr > sr_min) && (sr < sr_max)){
                     unit->sr = sr;
                     return i + 1;
              }

       }

       return 0;

}

/*
Serial worker: generates up to num_asu asymmetric units and writes them
as extended-XYZ blocks to out_path (truncated first). Returns the number
of asymmetric units actually generated.
*/
static int generate_asu_to_file(
    const char *out_path,
    molecule *mol,
    float sr_min,
    float sr_max,
    int *stoic,
    int n_mol_types,
    int num_asu,
    long max_attempts)
{
    int *n_atoms_in_mol = (int *)malloc(n_mol_types * sizeof(int));
    for(int i = 0; i < n_mol_types; i++)
        n_atoms_in_mol[i] = mol[i].num_of_atoms;

    asu unit;
    asu_init(&unit, stoic, n_atoms_in_mol, n_mol_types);

    FILE *out = fopen(out_path, "w");
    if(!out)
    {
        printf("***ERROR: cannot create %s\n", out_path);
        asu_free(&unit);
        free(n_atoms_in_mol);
        return 0;
    }

    long remain_attempts = max_attempts;
    int num_generated = 0;

    while(remain_attempts > 0 && num_generated < num_asu)
    {
        int budget = remain_attempts > INT_MAX ? INT_MAX : (int)remain_attempts;
        int used = generate_one_asu(&unit, mol, sr_min, sr_max, budget);
        if(used == 0)   // exhausted attempts without a valid placement
            break;
        remain_attempts -= used;
        num_generated++;
        write_asu_extxyz(&unit, out, num_generated);
    }

    fclose(out);
    asu_free(&unit);
    free(n_atoms_in_mol);
    return num_generated;
}

/*
Concatenates the per-rank shard files "<output_file>.rank<r>" into
output_file (in rank order) and deletes each shard afterwards.
*/
static void merge_shard_files(const char *output_file, int total_ranks)
{
    FILE *out = fopen(output_file, "w");
    if(!out)
    {
        printf("***ERROR: cannot create %s\n", output_file);
        return;
    }

    char shard_path[600];
    char buffer[8192];
    for(int r = 0; r < total_ranks; r++)
    {
        snprintf(shard_path, sizeof(shard_path), "%s.rank%d", output_file, r);
        FILE *shard = fopen(shard_path, "r");
        if(!shard)
            continue;
        size_t n;
        while((n = fread(buffer, 1, sizeof(buffer), shard)) > 0)
            fwrite(buffer, 1, n, out);
        fclose(shard);
        remove(shard_path);
    }
    fclose(out);
}

int generate_asymmetric_units(
    double *positions,
    int n_atoms_total,
    int ncols,
    char *species,
    int *n_atoms_per_mol,
    int n_mol_types,
    int *stoic,
    int n_mol_types_b,
    int num_structures,
    double sr_min,
    double sr_max,
    long max_attempts,
    int random_seed,
    char *output_file,
    MPI_Comm world_comm)
{
    int total_ranks, my_rank;
    MPI_Comm_size(world_comm, &total_ranks);
    MPI_Comm_rank(world_comm, &my_rank);

    if(ncols != 3)
    {
        if(my_rank == 0)
            printf("***ERROR: positions must have 3 columns, got %d\n", ncols);
        return 0;
    }
    if(n_mol_types != n_mol_types_b)
    {
        if(my_rank == 0)
            printf("***ERROR: n_atoms_per_mol and stoic must have equal length\n");
        return 0;
    }

    // Rebuild the per-type molecule array from the concatenated flat arrays.
    molecule *mol = (molecule *)malloc(n_mol_types * sizeof(molecule));
    int offset = 0;
    for(int m = 0; m < n_mol_types; m++)
    {
        int na = n_atoms_per_mol[m];
        mol[m].num_of_atoms = na;
        mol[m].X = (float *)malloc(na * sizeof(float));
        mol[m].Y = (float *)malloc(na * sizeof(float));
        mol[m].Z = (float *)malloc(na * sizeof(float));
        mol[m].atoms = (char *)malloc(2 * na * sizeof(char));
        for(int a = 0; a < na; a++)
        {
            int g = offset + a;
            mol[m].X[a] = (float)positions[g * 3 + 0];
            mol[m].Y[a] = (float)positions[g * 3 + 1];
            mol[m].Z[a] = (float)positions[g * 3 + 2];
            mol[m].atoms[2 * a]     = species[2 * g];
            mol[m].atoms[2 * a + 1] = species[2 * g + 1];
        }
        offset += na;
        recenter_molecule(&mol[m]);
    }

    // Seed the Mersenne-Twister per rank so each draws an independent stream.
    init_genrand((unsigned int)(random_seed + my_rank));

    // Split the requested count across ranks; low ranks absorb the remainder.
    int base = num_structures / total_ranks;
    int rem  = num_structures % total_ranks;
    int my_share = base + (my_rank < rem ? 1 : 0);

    // Each rank writes its structures to a private shard file.
    char shard_path[600];
    snprintf(shard_path, sizeof(shard_path), "%s.rank%d", output_file, my_rank);

    int my_generated = 0;
    if(my_share > 0)
        my_generated = generate_asu_to_file(shard_path, mol, (float)sr_min,
                                            (float)sr_max, stoic, n_mol_types,
                                            my_share, max_attempts);

    MPI_Barrier(world_comm);

    // Rank 0 merges the shards into the final output file.
    if(my_rank == 0)
        merge_shard_files(output_file, total_ranks);

    MPI_Barrier(world_comm);

    for(int m = 0; m < n_mol_types; m++)
    {
        free(mol[m].atoms);
        free(mol[m].X);
        free(mol[m].Y);
        free(mol[m].Z);
    }
    free(mol);

    int total_generated = 0;
    MPI_Allreduce(&my_generated, &total_generated, 1, MPI_INT, MPI_SUM, world_comm);

    if(my_rank == 0)
        printf("Generated %d/%d asymmetric units. Written to %s\n",
               total_generated, num_structures, output_file);

    return total_generated;
}
