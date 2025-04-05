"""Module of Shamirs Sharing and Leakage Resilient Sharing Schemes."""
import math
import numpy as np
import galois

# Large prime number to calculate the needed source length (Data Type: Long)
PRIME_MODULUS = 170141183460469231731687303715884105727

# Float using PRIME_MODULUS to calculate needed source length (Data Type: Float)
PRIME_MODULUS_BIT_LENGTH = math.ceil(math.log2(PRIME_MODULUS)) # (Value: ~127.0)



class ShamirSS:
    ''' Provides methods for encoding/decoding secret shares with Shamir's 
        secret sharing in a finite field. For education purposes only.
        Efficiency not considered!
    '''

    def __init__(self, p: int, w: int, rep: str="int"):
        ''' Initializes Galois Field in which operations take place 

            Arguments: 
            p -- prime characteristic of the field, typically p=2 
            w -- extension degree of field, so that the field has order 2^w 
            rep -- indicates how field elements are represented {"int", "power", "poly"}
        '''
        self.gf = galois.GF(p**w)
        self.gf.repr(rep)

    def deconstruct(self,
                    secret: int,
                    threshold: int,
                    num_shares: int,
                    rand_coeff: int=None) -> galois.FieldArray:
        ''' Split a secret into secret shares using Shamir's Secret Sharing scheme.

            Arguments: 
            secret -- Secret to be split. Must be a valid element of self.gf 
            threshold -- Reconstruction threshold (1 greater than polynomial degree) 
            num_shares -- Number of shares. Must have n > t 
            rand_coeff -- "random" coefficients to be used for polynomial. 
            If set to None, this method will randomly sample them.
        '''

        # TODO: better error handling, also assuming rand_coeff elements are in self.gf
        if rand_coeff is None:
            rand_coeff = self.gf.Random(threshold - 1)
        else:
            assert len(rand_coeff) == threshold - 1

        # NOTE: will raise ValueError if s does not belong to GF
        msg_array = np.append(self.gf(secret), rand_coeff)

        # Create a vandermonde matrix of valid elements in self.gf
        vandermonde = self.gf.Vandermonde(self.gf.primitive_element,
                                          num_shares,
                                          threshold)

        # Matrix Multiplication
        product = vandermonde @ msg_array

        shares = self.gf(np.zeros((num_shares, 2), dtype=int))
        shares[:,0] = vandermonde[:,1] # evaluation points
        shares[:,1] = product # value vector

        return shares

    def reconstruct(self,
                    shares: galois.FieldArray,
                    threshold: int=None,
                    num_shares: int=None) -> int:
        ''' Attempts to recover secret from given shares 
        
            Arguments: 
            shares -- two column array. First column is evaluation point,
                        second column is corresponding share value.
            threshold -- t value (may not be same value used in encoding). Assume t == n
            num_shares -- number of shares to use for decoding. n must be <= len(shares)

            NOTE: The sensible way of using this method, assuming t sufficient valid shares,
            is self.reconstruct(shares). If we have more than t sufficient valid shares,
            we can also do self.reconstruct(shares, t=t)
        '''

        # parameter validation and setup
        if num_shares is None or num_shares > len(shares) or num_shares < 1:
            num_shares = len(shares) # by default, try using all shares we are passed
        if threshold is None or threshold < 1:
            threshold = num_shares # by default, assume threshold is number of shares
        if threshold < num_shares:
            num_shares = threshold # we only need t shares to interpolate

        # NOTE: possible to have n < t, so our matrix will not have full rank.
        # Simulating: adversary aware they have insufficient shares and knows threshold
        # Nullspace: set of all possible solutions to compare to original message space

        # get first n shares and extract evaluation points / value vectors
        shares = shares[:num_shares]
        eval_points = shares[:,0]
        shares = shares[:,1]

        # build matrix for Ax = b
        matrix = self.gf(np.zeros((num_shares, threshold), dtype=int))
        for i, point in enumerate(eval_points):
            matrix[i] = self.gf([self.gf(point) ** j for j in range(0, threshold)])

        # TODO: Method not applicable to non-square matrices, try LU decomposition?
        x = np.linalg.solve(matrix, self.gf(shares))
        return x[0]


class LRSS:
    ''' Provides methods of encoding/decoding secret shares with Leakage Resilience. 
        Uses finite field elements.
    '''

    def __init__(self, p: int, w: int, rep: str="int"):
        ''' Initializes Galois Field in which operations take place 

            Arguments: 
            p -- prime characteristic of the field, typically p=2 
            w -- extension degree of field, so that the field has order 2^w 
            rep -- indicates how field elements are represented {"int", "power", "poly"}
        '''
        self.gf = galois.GF(p**w)
        self.gf.repr(rep)
        self.shamirs = ShamirSS(p, w, rep="int") # Same finite field as LRSS

    def inv_ext(self,
                share: galois.FieldArray,
                seed: galois.FieldArray,
                seed_length: int) -> galois.FieldArray:
        ''' Create a random source array from shares

            Arguments: 
            self -- Galois field inherited from LRSS Class 
            share -- two column array. First column is evaluation point,
                        second column is corresponding share value.
            seed -- "randomly" generated array of field elements. 
            seed_length -- length of seed array.
        '''

        source_array = self.gf.Random(len(seed))
        pivot = 0

        # Set pivot to index of first non-zero seed value.
        for i in range(seed_length):
            if seed[i] != self.gf(0):
                pivot = i # Prevents division by zero
                break

        share_value = self.gf(share[1])

        # Let y equal the value of the share.
        # Let w equal source value at its index (w_0 is the first value of the source)
        # Let s equal seed value at its index (s_0 is the first value of the source)
        #
        # Such that: y = w_0 * s_0 + w_1 * s_1 ... w_d * s_d, where d equals seed length
        #
        # Solving for w_0: w_0 = ((y - (w_1 * s_1)) // s_0), for d = 1.
        #
        # To perform the subtraction term: (y - (w_1 * s_1)),
        # For every index of the seed that is not the pivot, perform subtraction.
        # The value will decrement by the product of source_array and seed's values at each index.
        for i in range(seed_length):
            if i != pivot:
                share_value = share_value - (source_array[i]*seed[i])

        # Create source at pivot index.
        source_array[pivot] = share_value/(seed[pivot])

        return source_array

    def ext(self, source: int, seed: galois.FieldArray) -> int:
        ''' Calculates a dot product from a seed and source.

            Arguments: 
            self -- Galois field inherited from LRSS Class 
            source -- semirandomly generated field element array based on share values.
            seed -- "randomly" generated array of field elements. 
        '''

        dot_product = np.dot(self.gf(source), self.gf(seed))

        return dot_product


    def lr_share(self,
                 secret: int,
                 threshold: int,
                 num_shares: int,
                 lbudget: float,
                 lerror: float) -> galois.FieldArray:
        ''' Split a secret into secret shares using a leakage resilient sharing scheme.

            Arguments: 
            self - The designated Galois Field
            secret -- Secret to be split. Must be a valid element of self.gf 
            threshold -- Reconstruction threshold (1 greater than polynomial degree)
            num_shares -- Number of shares. Must have n > t 
            lbudget -- The leakage of bits allowed to be acquired by the adversary.
            lerror -- A parameter accounting for potential error due to leakage.
        '''

        # Find source length.
        ext_error = lerror / (6 * num_shares) 
        min_entropy = PRIME_MODULUS_BIT_LENGTH + 2 * np.log2(1 / ext_error) - 2
        source_length = np.ceil((min_entropy + lbudget) / PRIME_MODULUS_BIT_LENGTH)

        seed_length = int(source_length)

        # Perform Shamir's Secret Sharing
        msg_shares = self.shamirs.deconstruct(secret, threshold, num_shares, None)

        seed = self.gf.Random(seed_length)
        while seed.all() == self.gf(0):
            seed = self.gf.Random(seed_length) # Ensure non-zero seed

        # Secret Share seed elements
        seed_shares = self.gf(np.zeros((seed_length, num_shares, 2),dtype=int))
        for i in range(seed_length):
            seed_shares[i] = self.shamirs.deconstruct(seed[i], 2, num_shares, None)

        source_array = []
        share = []

        # Create share through sources and seed shares
        for i in range(num_shares):
            source_array.insert(i, self.inv_ext(msg_shares[i], seed, seed_length))
            share.insert(i, (source_array[i], seed_shares[:,i]))

        return share

    def lr_share_rec(self,
                     shares: galois.FieldArray,
                     threshold: int,
                     num_shares: int) -> int:
        ''' Reconstruct shares into secret using a leakage resilient sharing scheme.

            Arguments: 
            self - The designated Galois Field
            shares -- Shares to be reconstructed. Must contain valid elements of self.gf 
            threshold -- Reconstruction threshold (1 greater than polynomial degree) 
            num_shares -- Number of shares. Must have n > t 
            metadata -- Currently NONE, but could hold additional information.
        '''

        seed_length = len(shares[0][0]) # number of seed shares
        seed_vec = self.gf(np.zeros((seed_length, num_shares , 2), dtype=int))
        source_list = []
        extractor_output = []

        for i in range(seed_length):
            for j, share in enumerate(shares):
                seed_vec[i][j] = share[1][i] # seed shares from evaluation points

        seed = self.gf.Zeros(seed_length)

        # NOTE: The threshold for seed reconstruction is 2, as performed in the paper,
        # "Short Leakage Resilient and Non-malleable Secret Sharing Schemes".
        for i in range(seed_length):
            seed[i] = self.shamirs.reconstruct(seed_vec[i], 2, num_shares)

        eval_points = []

        for i in range(threshold):
            eval_points.insert(i, shares[i][1][0][0])
            source_list.insert(i, shares[i][0])
            extractor_output.insert(i, self.ext(source_list[i], seed))

        # NOTE: vector is transposed to fit reconstruction algorithm
        rec_vector = np.transpose(np.array((self.gf(eval_points),
                                            self.gf(extractor_output))))

        msg = self.shamirs.reconstruct(rec_vector, threshold, num_shares)

        return msg
