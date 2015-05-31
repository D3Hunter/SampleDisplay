#ifndef _RANDOM_H_
#define _RANDOM_H_

// maximum allowed space dimension
#define SOBOL_MAX_DIMENSION 40

// bit count; assumes sizeof(int) >= 32-bit
#define SOBOL_BIT_COUNT 30

// degrees of the primitive polynomials
const int degree_table[SOBOL_MAX_DIMENSION] = {
	0, 1, 2, 3, 3, 4, 4, 5, 5, 5,
	5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
	7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
	7, 7, 7, 7, 7, 7, 7, 8, 8, 8
};

// primitive polynomials in binary encoding
const int primitive_polynomials[SOBOL_MAX_DIMENSION] =	{
	1,     3,   7,  11,  13,  19,  25,  37,  59,  47,
	61,   55,  41,  67,  97,  91, 109, 103, 115, 131,
	193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
	213, 191, 253, 203, 211, 239, 247, 285, 369, 299
};
// initial values for direction tables, following
// Bratley+Fox, taken from [Sobol+Levitan, preprint 1976]
const int v_init[8][SOBOL_MAX_DIMENSION] = {
	{
		0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1
	},
	{
		0, 0, 1, 3, 1, 3, 1, 3, 3, 1,
			3, 1, 3, 1, 3, 1, 1, 3, 1, 3,
			1, 3, 1, 3, 3, 1, 3, 1, 3, 1,
			3, 1, 1, 3, 1, 3, 1, 3, 1, 3
		}, 
		{
			0, 0, 0, 7, 5, 1, 3, 3, 7, 5,
				5, 7, 7, 1, 3, 3, 7, 5, 1, 1,
				5, 3, 3, 1, 7, 5, 1, 3, 3, 7,
				5, 1, 1, 5, 7, 7, 5, 1, 3, 3
		}, 
		{
			0,  0,  0,  0,  0,  1,  7,  9, 13, 11,
				1,  3,  7,  9,  5, 13, 13, 11,  3, 15,
				5,  3, 15,  7,  9, 13,  9,  1, 11,  7,
				5, 15,  1, 15, 11,  5,  3,  1,  7,  9
			}, 
			{
				0,  0,  0,  0,  0,  0,  0,  9,  3, 27,
					15, 29, 21, 23, 19, 11, 25,  7, 13, 17,
					1, 25, 29,  3, 31, 11,  5, 23, 27, 19,
					21,  5,  1, 17, 13,  7, 15,  9, 31,  9
			}, 
			{
				0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
					0,  0,  0, 37, 33,  7,  5, 11, 39, 63,
					27, 17, 15, 23, 29,  3, 21, 13, 31, 25,
					9, 49, 33, 19, 29, 11, 19, 27, 15, 25
				}, 
				{
					0,   0,  0,  0,  0,  0,    0,  0,  0,   0,
						0,   0,  0,  0,  0,  0,    0,  0,  0,  13,
						33, 115, 41, 79, 17,  29, 119, 75, 73, 105,
						7,  59, 65, 21,  3, 113,  61, 89, 45, 107
				}, 
				{
					0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
						0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
						0, 0, 0, 0, 0, 0, 0, 0,  0,  0,
						0, 0, 0, 0, 0, 0, 0, 7, 23, 39
					}
};

///////////////////////////////////////////////////////////////////////
// Class				:	CSobol
// Description			:	Sobol quasi random generator
// Comments				:
template <int dimension> class CSobol {
public:
	CSobol(int seq = 1) {
		init(seq);
	}

	void	init(int seq = 1) {
		unsigned int	i_dim;
		int				j, k;
		int				ell;

		//assert(dimension < SOBOL_MAX_DIMENSION);

		// Initialize direction table in dimension 0.
		for(k=0; k<SOBOL_BIT_COUNT; k++) v_direction[k][0] = 1;

		// Initialize in remaining dimensions.
		for(i_dim=1; i_dim<dimension; i_dim++) {
			const int	poly_index	= i_dim;
			const int	degree_i	= degree_table[poly_index];
			int			includ[8];

			// Expand the polynomial bit pattern to separate
			// components of the logical array includ[].
			int p_i = primitive_polynomials[poly_index];
			for(k = degree_i-1; k >= 0; k--) {
				includ[k] = ((p_i % 2) == 1);
				p_i >>= 1;
			}

			// Leading elements for dimension i come from v_init[][].
			for(j=0; j<degree_i; j++) v_direction[j][i_dim] = v_init[j][i_dim];

			// Calculate remaining elements for this dimension,
			// as explained in Bratley+Fox, section 2.
			for(j=degree_i; j<SOBOL_BIT_COUNT; j++) {
				int newv = v_direction[j-degree_i][i_dim];
				ell = 1;
				for(k=0; k<degree_i; k++) {
					ell *= 2;
					if(includ[k]) newv ^= (ell * v_direction[j-k-1][i_dim]);
				}
				v_direction[j][i_dim] = newv;
			}
		}

		// Multiply columns of v by appropriate power of 2.
		ell = 1;
		for(j=SOBOL_BIT_COUNT-1-1; j>=0; j--) {
			ell *= 2;
			for(i_dim=0; i_dim<dimension; i_dim++) {
				v_direction[j][i_dim] *= ell;
			}
		}

		// 1/(common denominator of the elements in v_direction)
		last_denominator_inv = (float) (1.0 /(2.0 * ell));

		// final setup 
		sequence_count = seq;
		for(i_dim=0; i_dim<dimension; i_dim++) last_numerator_vec[i_dim] = 0;
	}



	void			get(float *v) {
		unsigned int i_dimension;

		// Find the position of the least-significant zero in count.
		int ell = 0;
		int c = sequence_count;
		while(1) {
			++ell;
			if((c % 2) == 1) c >>= 1;
			else break;
		}

		for(i_dimension=0; i_dimension<dimension; i_dimension++) {
			const int direction_i     = v_direction[ell-1][i_dimension];
			const int old_numerator_i = last_numerator_vec[i_dimension];
			const int new_numerator_i = old_numerator_i ^ direction_i;
			last_numerator_vec[i_dimension] = new_numerator_i;
			v[i_dimension] = new_numerator_i * last_denominator_inv;
		}

		sequence_count++;
		if (sequence_count >= (1 << SOBOL_BIT_COUNT)) {
			sequence_count	=	0;
		}
	}

	unsigned int	sequence_count;
	float			last_denominator_inv;
	int				last_numerator_vec[SOBOL_MAX_DIMENSION];
	int				v_direction[SOBOL_BIT_COUNT][SOBOL_MAX_DIMENSION];
};

#endif
