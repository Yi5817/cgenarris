#include "cgenarris/asu_generation.h"
#include "asu_utils.h"
#include "cgenarris/asu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "algebra.h"
#include "cgenarris/molecule.h"
#include "molecule_utils.h"
#include "randomgen.h"
#include "check_structure.h"   // find_mol_len

#define ASU_PATH_MAX 4096

float asu_box_length(const molecule *mol, const int *stoic, int n_mol_types)
{
    float volume = 0;
    for(int m = 0; m < n_mol_types; m++)
    {
        float d = find_mol_len(mol[m].X, mol[m].Y, mol[m].Z, mol[m].num_of_atoms);
        volume += d * d * d * stoic[m];
    }
    return cbrtf(3.0f * volume);
}

void asu_place_random(asu *unit, const molecule *mol, float box_len)
{
    float rot[3][3];
    float trans[3];

    for(int k = 0; k < unit->n_mols; k++)
    {
        const molecule *m = mol + unit->mol_types[k];
        int first = unit->mol_index[k];

        generate_random_rotation_matrix(rot);
        generate_random_translation_vector(trans);   // fractional, in [0, 1)

        for(int a = 0; a < m->num_of_atoms; a++)
        {
            float atom[3] = {m->X[a], m->Y[a], m->Z[a]};
            vector3_mat3b3_multiply(rot, atom, atom);
            unit->Xcord[first + a] = atom[0] + box_len * trans[0];
            unit->Ycord[first + a] = atom[1] + box_len * trans[1];
            unit->Zcord[first + a] = atom[2] + box_len * trans[2];
        }
    }
}

float asu_min_sr(const asu *unit, float stop_below)
{
    float min_ratio_sq = HUGE_VALF;   // (d_ij / (r_i + r_j))^2; n_mols >= 2
    float stop_sq = stop_below * stop_below;

    for(int p = 0; p < unit->n_mols - 1; p++)
    {
        int p0 = unit->mol_index[p];
        int p1 = p0 + unit->n_atoms_in_mol[unit->mol_types[p]];
        for(int q = p + 1; q < unit->n_mols; q++)
        {
            int q0 = unit->mol_index[q];
            int q1 = q0 + unit->n_atoms_in_mol[unit->mol_types[q]];
            for(int i = p0; i < p1; i++)
            {
                for(int j = q0; j < q1; j++)
                {
                    float dx = unit->Xcord[i] - unit->Xcord[j];
                    float dy = unit->Ycord[i] - unit->Ycord[j];
                    float dz = unit->Zcord[i] - unit->Zcord[j];
                    float rsum = unit->vdw_radii[i] + unit->vdw_radii[j];
                    float ratio_sq = (dx * dx + dy * dy + dz * dz) / (rsum * rsum);
                    if(ratio_sq < min_ratio_sq)
                        min_ratio_sq = ratio_sq;
                    if(ratio_sq < stop_sq)
                        return sqrtf(ratio_sq);   // clash: the unit is rejected anyway
                }
            }
        }
    }
    return sqrtf(min_ratio_sq);
}

long asu_generate_one(asu *unit, const molecule *mol, float box_len,
                      float sr_min, float sr_max, long max_attempts)
{
    for(long attempt = 1; attempt <= max_attempts; attempt++)
    {
        asu_place_random(unit, mol, box_len);
        float sr = asu_min_sr(unit, sr_min);
        if(sr > sr_min && sr < sr_max)
        {
            unit->sr = sr;
            asu_recenter(unit);
            return attempt;
        }
    }
    return 0;
}

/*
Serial worker: generates up to num_asu asymmetric units and writes them as
geometry.out blocks numbered from first_number to out_path (truncated
first). Returns the number of units written (or -1 if the file cannot be
created) and adds the number of placements tried to *attempts.
*/
static int generate_to_file(const char *out_path, asu *unit, const molecule *mol,
                            float box_len, float sr_min, float sr_max,
                            int num_asu, int first_number, long max_attempts,
                            long *attempts)
{
    FILE *out = fopen(out_path, "w");
    if(!out)
    {
        fprintf(stderr, "***ERROR: asu_generate_mpi: cannot create %s\n", out_path);
        return -1;
    }

    int n_written = 0;
    while(n_written < num_asu)
    {
        long used = asu_generate_one(unit, mol, box_len, sr_min, sr_max, max_attempts);
        if(used == 0)
        {
            *attempts += max_attempts;
            fprintf(stderr, "**WARNING: no asymmetric unit found within "
                    "%ld attempts; stopping after %d\n", max_attempts, n_written);
            break;
        }
        *attempts += used;
        n_written++;
        asu_write_block(unit, out, first_number + n_written - 1);
    }

    fclose(out);
    return n_written;
}

/*
Concatenates the per-rank shard files "<output_file>.rank<r>" into
output_file (in rank order) and deletes each shard afterwards. Each shard
numbers its blocks from 1; the "#structure_number" lines are rewritten so
the merged file is numbered consecutively even when a rank stopped early.
*/
static int merge_shard_files(const char *output_file, int total_ranks)
{
    char shard_path[ASU_PATH_MAX];
    for(int r = 0; r < total_ranks; r++)
    {
        snprintf(shard_path, sizeof(shard_path), "%s.rank%d", output_file, r);
        FILE *shard = fopen(shard_path, "r");
        if(!shard)
        {
            fprintf(stderr, "***ERROR: asu_generate_mpi: cannot read %s; "
                    "output_file must be on a filesystem shared by all ranks\n",
                    shard_path);
            return -1;
        }
        fclose(shard);
    }

    FILE *out = fopen(output_file, "w");
    if(!out)
    {
        fprintf(stderr, "***ERROR: asu_generate_mpi: cannot create %s\n", output_file);
        return -1;
    }

    static const char number_tag[] = "#structure_number = ";
    char line[8192];
    int structure_number = 0;
    for(int r = 0; r < total_ranks; r++)
    {
        snprintf(shard_path, sizeof(shard_path), "%s.rank%d", output_file, r);
        FILE *shard = fopen(shard_path, "r");
        if(!shard)
            continue;   // checked above; only a concurrent removal gets here
        int at_line_start = 1;
        while(fgets(line, sizeof(line), shard))
        {
            if(at_line_start && !strncmp(line, number_tag, sizeof(number_tag) - 1))
                fprintf(out, "%s%d\n", number_tag, ++structure_number);
            else
                fputs(line, out);
            at_line_start = line[strlen(line) - 1] == '\n';
        }
        fclose(shard);
        remove(shard_path);
    }
    fclose(out);
    return 0;
}

// Validates the user-facing arguments; prints on rank 0 only.
static int check_arguments(int my_rank, int n_mol_types, int num_structures,
                           double sr_min, double sr_max, long max_attempts,
                           const char *output_file)
{
    const char *problem = NULL;
    if(n_mol_types < 1)
        problem = "n_mol_types must be >= 1";
    else if(num_structures < 0)
        problem = "num_structures must be >= 0";
    else if(!(sr_min >= 0 && sr_min < sr_max))
        problem = "sr window must satisfy 0 <= sr_min < sr_max";
    else if(max_attempts < 1)
        problem = "max_attempts must be >= 1";
    else if(!output_file || output_file[0] == '\0')
        problem = "output_file must not be empty";
    else if(strlen(output_file) + 16 > ASU_PATH_MAX)
        problem = "output_file path is too long";

    if(problem)
    {
        if(my_rank == 0)
            fprintf(stderr, "***ERROR: asu_generate_mpi: %s\n", problem);
        return -1;
    }
    return 0;
}

int asu_generate_mpi(molecule *mol, int n_mol_types, const int *stoichiometry,
                     int num_structures, double sr_min, double sr_max,
                     long max_attempts, int random_seed,
                     const char *output_file, MPI_Comm comm)
{
    int total_ranks, my_rank;
    MPI_Comm_size(comm, &total_ranks);
    MPI_Comm_rank(comm, &my_rank);

    if(check_arguments(my_rank, n_mol_types, num_structures, sr_min, sr_max,
                       max_attempts, output_file))
        return -1;

    asu unit;
    if(asu_init(&unit, mol, stoichiometry, n_mol_types))
        return -1;

    for(int m = 0; m < n_mol_types; m++)
        recenter_molecule(mol + m);
    float box_len = asu_box_length(mol, stoichiometry, n_mol_types);

    // Seed the Mersenne Twister per rank so each draws an independent stream.
    if(random_seed == 0)
    {
        if(my_rank == 0)
            random_seed = (int)(time(NULL) % 1000000000L) + 1;
        MPI_Bcast(&random_seed, 1, MPI_INT, 0, comm);
    }
    init_genrand((unsigned int)(random_seed + my_rank));

    if(my_rank == 0)
    {
        printf("ASYMMETRIC UNIT GENERATION:\n");
        printf("-----------------------------\n");
        printf("Molecule types:                %d\n", n_mol_types);
        printf("Stoichiometry:                 ");
        for(int m = 0; m < n_mol_types; m++)
            printf("%d%s", stoichiometry[m], m + 1 < n_mol_types ? ":" : "\n");
        printf("Atoms per asymmetric unit:     %d\n", unit.n_atoms);
        printf("Placement box edge (Angstrom): %.3f\n", box_len);
        printf("sr window:                     (%.3f, %.3f)\n", sr_min, sr_max);
        printf("Max attempts per unit:         %ld\n", max_attempts);
        printf("Random seed:                   %d\n", random_seed);
        printf("Requested units:               %d on %d ranks\n",
               num_structures, total_ranks);
        printf("-----------------------------\n");
        fflush(stdout);
    }

    // Split the requested count across ranks; low ranks absorb the remainder.
    // Each shard is numbered from 1; merge_shard_files renumbers consecutively.
    int base = num_structures / total_ranks;
    int rem = num_structures % total_ranks;
    int my_share = base + (my_rank < rem ? 1 : 0);

    char shard_path[ASU_PATH_MAX];
    snprintf(shard_path, sizeof(shard_path), "%s.rank%d", output_file, my_rank);

    long my_attempts = 0;
    int my_written = generate_to_file(shard_path, &unit, mol, box_len,
                                      (float)sr_min, (float)sr_max, my_share,
                                      1, max_attempts, &my_attempts);
    asu_free(&unit);

    int my_error = my_written < 0;
    int any_error = 0;
    MPI_Allreduce(&my_error, &any_error, 1, MPI_INT, MPI_MAX, comm);

    if(my_rank == 0 && !any_error)
        any_error = merge_shard_files(output_file, total_ranks) != 0;
    MPI_Bcast(&any_error, 1, MPI_INT, 0, comm);
    if(any_error)
    {
        remove(shard_path);   // leave no partial shards behind
        return -1;
    }

    int total_written = 0;
    long total_attempts = 0;
    MPI_Allreduce(&my_written, &total_written, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&my_attempts, &total_attempts, 1, MPI_LONG, MPI_SUM, comm);

    if(my_rank == 0)
    {
        printf("Generated %d/%d asymmetric units in %ld attempts "
               "(acceptance %.4f). Written to %s\n",
               total_written, num_structures, total_attempts,
               total_attempts ? (double)total_written / total_attempts : 0.0,
               output_file);
        fflush(stdout);
    }
    return total_written;
}

int generate_asymmetric_units(
    double *positions,
    int n_atoms_total,
    int ncols,
    char *species,
    int *n_atoms_per_mol,
    int n_mol_types,
    int *stoichiometry,
    int n_stoichiometry,
    int num_structures,
    double sr_min,
    double sr_max,
    long max_attempts,
    int random_seed,
    char *output_file,
    MPI_Comm world_comm)
{
    int my_rank;
    MPI_Comm_rank(world_comm, &my_rank);

    int atom_sum = 0;
    int min_atoms = 1;
    for(int m = 0; m < n_mol_types; m++)
    {
        atom_sum += n_atoms_per_mol[m];
        if(n_atoms_per_mol[m] < min_atoms)
            min_atoms = n_atoms_per_mol[m];
    }

    const char *problem = NULL;
    if(ncols != 3)
        problem = "positions must have 3 columns";
    else if(n_mol_types < 1)
        problem = "n_atoms_per_mol must not be empty";
    else if(min_atoms < 1)
        problem = "every entry of n_atoms_per_mol must be >= 1";
    else if(n_stoichiometry != n_mol_types)
        problem = "stoichiometry and n_atoms_per_mol must have equal length";
    else if(atom_sum != n_atoms_total)
        problem = "sum(n_atoms_per_mol) must equal the number of positions";
    else if(!species || strlen(species) != 2 * (size_t)n_atoms_total)
        problem = "species must hold exactly 2 chars per atom";

    if(problem)
    {
        if(my_rank == 0)
            fprintf(stderr, "***ERROR: generate_asymmetric_units: %s\n", problem);
        return -1;
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
    }

    int total = asu_generate_mpi(mol, n_mol_types, stoichiometry, num_structures,
                                 sr_min, sr_max, max_attempts, random_seed,
                                 output_file, world_comm);

    for(int m = 0; m < n_mol_types; m++)
    {
        free(mol[m].atoms);
        free(mol[m].X);
        free(mol[m].Y);
        free(mol[m].Z);
    }
    free(mol);
    return total;
}
