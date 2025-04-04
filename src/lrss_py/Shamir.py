# Required imports for Shamir Module:
import numpy as np
import galois


# Global Variables:

# PRIME_MODULUS is a large prime number used to calculate the needed source length (Data Type: Float) 
PRIME_MODULUS = float(170141183460469231731687303715884105727)

# PRIME_MODULUS_BIT_LENGTH is a float that uses PRIME_MODULUS used to calculate the needed source length (Data Type: Float) 
# (Value: 127.0)
PRIME_MODULUS_BIT_LENGTH = np.ceil(np.log2(PRIME_MODULUS))



class ShamirSS:
    ''' Provides methods for encoding/decoding secret shares with Shamir's 
        secret sharing in a finite field. For education purposes only.
        Efficiency not considered!
    '''


    def __init__(self, p, w, repr="int"):
        ''' Initializes Galois Field in which operations take place 

            Arguments: 
            p -- prime characteristic of the field, typically p=2 
            w -- extension degree of field, so that the field has order 2^w 
            repr -- one of {"int", "poly", "power"}, indicating how elements
                    of the field should be represented.
        '''
        self.GF = galois.GF(p**w)
        self.GF.repr(repr)

    def deconstruct(self, secret, threshold, number_of_shares, rand_coeff=None):
        ''' Split a secret into secret shares using Shamir's Secret Sharing scheme.

            Arguments: 
            secret -- Secret to be split. Must be a valid element of self.GF 
            threshold -- Reconstruction threshold (1 greater than polynomial degree) 
            number_of_shares -- Number of shares. Must have n > t 
            rand_coeff -- "random" coefficients to be used for polynomial. 
            If set to None, this method will randomly sample them.
        ''' 

        # Check rand_coeff parameter. If None, then assign it as a random element of self.GF
        if (rand_coeff == None):
            rand_coeff = self.GF.Random(threshold - 1)
        # Else, ensure that the length of rand_coeff is one less than the threshold value.
        else:
            assert len(rand_coeff) == threshold - 1 # TODO: better error handling, also assuming rand_coeff elements are in self.GF

        # Append the rand_coeff value to the secret in an array (msg). 
        msg = np.append(self.GF(secret), rand_coeff) # NOTE: will raise ValueError if s does not belong to GF

        # Create a vandermonde matrix of valid elements in self.GF
        vandermonde = self.GF.Vandermonde(self.GF.primitive_element,
                                          number_of_shares,
                                          threshold)

        # Perform matrix multiplication on the vandermonde matrix and the msg array
        product = vandermonde @ msg

        # Set a new variable, shares, to be a zero matrix of size (n,2)
        shares = self.GF(np.zeros((number_of_shares, 2), dtype=int))

        # Set shares' first column to have the first column of elements from the vandermonde (will serve as the evaluation points)
        shares[:,0] = vandermonde[:,1] # eval points

        # Set shares'second column to be the product of the matrix multiplication (will serve as the value vector)
        shares[:,1] = product

        # Output the new shares array
        return shares
    
    def reconstruct(self, shares, t=None, n=None):
        ''' Attempts to recover secret from given shares 
        
            Arguments: 
            shares -- two column array. First column is evaluation point,
                        second column is corresponding share value.
            t -- threshold value (may not be same value used in encoding). By default assume t == n
            n -- number of shares to use for decoding. n must be <= len(shares)

            NOTE: The sensible way of using this method, assuming t sufficient valid shares,
            is self.reconstruct(shares). If we have more than t sufficient valid shares, we can also do
            self.reconstruct(shares, t=t)
        '''
        
        # parameter validation and setup
        if (n == None or n > len(shares) or n < 1):
            n = len(shares) # by default, try using all shares we are passed
        if (t == None or t < 1):
            t = n # by default, assume threshold is number of shares
        if (t < n):
            n = t # if t < n, we only need t shares to interpolate, so set n = t

        # NOTE: at this point it is possible to have n < t, so our matrix will not have full rank.
        # This case makes sense for simulating when adversary is aware they have insufficient shares and knows true threshold.
        # Can look at null space for set of all possible solutions (and then compare this to original message space)
        
        # get first n shares and extract evaluation points / values
        shares = shares[:n]
        eval_points = shares[:,0]
        shares = shares[:,1]


        # build n by t coefficient matrix using eval_points
        A = self.GF(np.zeros((n,t), dtype=int))
        for i, point in enumerate(eval_points):
            A[i] = self.GF([self.GF(point) ** j for j in range(0, t)])

        # solve Ax = b (TODO: np.linalg.solve does not work for non-square matrices, maybe use LU decomposition?
        x = np.linalg.solve(A, self.GF(shares))

        # return first element of x (the secret)
        return x[0]


class LRSS:
    ''' Provides methods for encoding/decoding secret shares with Leakage Resilience Secret Sharing in a finite field. 
    '''

    def __init__(self, p, w, repr="int"):
        ''' Initializes Galois Field in which operations take place 

            Arguments: 
            p -- prime characteristic of the field, typically p=2 
            w -- extension degree of field, so that the field has order 2^w 
            repr -- one of {"int", "poly", "power"}, indicating how elements
                    of the field should be represented.
        '''
        self.GF = galois.GF(p**w)
        self.GF.repr(repr)
        self.Shamirs = ShamirSS(p,w)
    
    def InvExt(self, share, threshold, number_of_shares, seed, seed_length):
        ''' Create a random source array from shares

            Arguments: 
            self -- Galois field inherited from LRSS Class 
            share -- two column array. First column is evaluation point,
                        second column is corresponding share value.
            threshold -- threshold value (may not be same value used in encoding). By default assume t == n.
            number_of_shares -- Total number of shares. Must have n > t 
            seed -- "randomly" generated array of field elements. 
            seed_length -- length of seed array.
        ''' 
        # parameter validation and setup
        if (number_of_shares == None or number_of_shares > len(share) or number_of_shares < 1):
            number_of_shares = len(share) # by default, try using all shares we are passed
        if (threshold == None or threshold < 1):
            threshold = number_of_shares # by default, assume threshold is number of shares
        if (threshold < number_of_shares):
            number_of_shares = threshold # if t < n, we only need t shares to interpolate, so set n = t

        # NOTE: at this point it is possible to have n < t, so our matrix will not have full rank.
        # This case makes sense for simulating when adversary is aware they have insufficient shares and knows true threshold.
        # Can look at null space for set of all possible solutions (and then compare this to original message space)
        
        # Create an array for the source elements that is the length of the seed and fill it with random field elements.
        source_array = self.GF.Random(len(seed))  

        # Create a variable for the pivot. The pivot identifies the first non-zero element in the seed (Prevents undefined quotient).
        pivot = 0

        # Iterate through the seed, looking for the first non-zero value. Set the pivot to be that value's index.
        for i in range(seed_length):
            if seed[i] != self.GF(0):
                pivot = i
                break

        # Create a variable for the share's value.
        share_value = self.GF(share[1])


        # Let y equal the value of the share.
        # Let w be the source value at its respective index (w_0 is the first value of the source)
        # Let s be the seed value at its respective index (s_0 is the first value of the source)
        #
        # Such that: y = w_0 * s_0 + w_1 * s_1 ... w_d * s_d, where d is the length of the seed
        #
        # Solving for w_0: w_0 = ((y - (w_1 * s_1)) // s_0), for d = 1.
        #
        # To perform the subtraction term: (y - (w_1 * s_1)),
        # For every index of the seed that is not the pivot index, decrement the value of share_value. 
        # The value will decrement by the product of source_array and seed's values at each index.

        for i in range(seed_length):
            if i != pivot:
                share_value = share_value - (source_array[i]*seed[i])

        # Now we divide that term by the seed at the pivot index, which solves for the source at the pivot index.
        source_array[pivot] = share_value/(seed[pivot])

        # We then return the entire source array (OUTPUT of InvExt).
        return source_array

    def ext(self, source, seed):
        ''' Calculates a dot product from a seed and source.

            Arguments: 
            self -- Galois field inherited from LRSS Class 
            source -- semirandomly generated field element array based on share values.
            seed -- "randomly" generated array of field elements. 
        ''' 
        # Take the dot product of the source and seed values.
        dot_product = np.dot(self.GF(source), self.GF(seed))

        # Return the dot product (OUTPUT of ext)
        return dot_product


    def lr_share(self, secret, threshold, number_of_shares, lbudget, lerror):
        ''' Split a secret into secret shares using a leakage resilient sharing scheme.

            Arguments: 
            self - The designated Galois Field
            secret -- Secret to be split. Must be a valid element of self.GF 
            threshold -- Reconstruction threshold (1 greater than polynomial degree)
            number_of_shares -- Number of shares. Must have n > t 
            lbudget -- The leakage of bits allowed to be acquired by the adversary.
            lerror -- A parameter accounting for potential error due to leakage.
        ''' 

        # Calculate the source length.
        ext_error = lerror / (6 * number_of_shares) 
        min_entropy = PRIME_MODULUS_BIT_LENGTH + 2 * np.log2(1 / ext_error) - 2
        source_length = np.ceil((min_entropy + lbudget) / PRIME_MODULUS_BIT_LENGTH)

        # Set seed length to the source length.
        seed_length = int(source_length)

        # Deconstruct the secret into n message shares using the function: deconstruct.
        msg_shares = self.Shamirs.deconstruct(secret, threshold, number_of_shares, None)

        # Create a seed that is an array of random elements within the galois field and of the same length as the source length.
        seed = self.GF.Random(seed_length)
        
        # Ensure that atleast one element of the seed is nonzero, if not, recreate seed array.
        while seed.all() == self.GF(0):
            seed = self.GF.Random(seed_length)

        # Create a zero matrix to store seed shares of size (seed_length, number of shares, 2)
        seed_shares = self.GF(np.zeros((seed_length, number_of_shares, 2),dtype=int))

        # Enumerate through seed length, for every index in the seed array, perform deconstruction on that index's element.
        # Assign the deconstruction result to the same index in the seed matrix.
        for i in range(seed_length):
            seed_shares[i] = self.Shamirs.deconstruct(seed[i], 2, number_of_shares, None)


        # Create an empty list for the source array
        source_array = []

        # Create an empty list for the final share array
        share = []

        # Index through each share via for loop.
        for i in range(number_of_shares):

            # Fill each index of the source array with the output from the Inverse Extractor at each individual message share.
            source_array.insert(i,self.InvExt(msg_shares[i], threshold, number_of_shares, seed, seed_length))

            # Fill the final share array with a tuple containting the source at index i and the value vector of seed share array.
            share.insert(i, (source_array[i], seed_shares[:,i]))

        # Return the final share array (OUTPUT of lr_share).
        return share

    def lr_share_rec(self, shares, threshold, number_of_shares, metadata):
        ''' Reconstruct secret shares into the original secret using a leakage resilient sharing scheme.

            Arguments: 
            self - The designated Galois Field
            shares -- Shares to be reconstructed. Must contain valid elements of self.GF 
            threshold -- Reconstruction threshold (1 greater than polynomial degree) 
            number_of_shares -- Number of shares. Must have n > t 
            metadata -- Currently NONE, but could hold additional information.
        ''' 

        # Set the seed length, which is equivalent to the number of the seed shares. 
        seed_length = len(shares[0][0])

        # Create a 3-dimnsional zero vector of size seed_length by number_of_shares by 2 called seed_vec.
        seed_vec = self.GF(np.zeros((seed_length, number_of_shares , 2), dtype=int))

        # Create an empty list for the source array.
        source_list = []

        # Create an empty list for the extractor output.
        extractor_output = []

        # Insert the seed shares from each set of equilibrium points into seed_vec.
        for i in range(seed_length):
            for j in range(len(shares)):
                seed_vec[i][j] = shares[j][1][i]

        # Create a variable array, seed, to hold the reconstructed seed (size of seed_length)
        seed = self.GF.Zeros(seed_length)

        # Perform reconstruction on each set of seed shares from seed_vec and insert them into seed.
        # The threshold for seed reconstruction is 2, as performed in the paper:
        ## "Short Leakage Resilient and Non-malleable Secret Sharing Schemes".
        for i in range(seed_length):
            seed[i] = self.Shamirs.reconstruct(seed_vec[i], 2, number_of_shares)

        # Create an empty list for the evaluation points.
        eval_points = []

        # Indexing up to the threshold amount:
        for i in range(threshold):

            # Insert the evaluation points that were used for the seed_shares into eval_points.
            eval_points.insert(i, shares[i][1][0][0])

            # Insert the source arrays into source_list.
            source_list.insert(i, shares[i][0])
            
            # Insert the outputs of the Extractor using source array with the reconstructed seed into extractor_output.
            extractor_output.insert(i, self.ext(source_list[i], seed))
        
        # Create the vector that will be used for reconstruction: a transposed array of the evaluation points with the extractor outputs.
        rec_vector = np.transpose(np.array((self.GF(eval_points), self.GF(extractor_output))))

        # Solve for the original message via reconstruction.
        msg = self.Shamirs.reconstruct(rec_vector, threshold, number_of_shares)

        # Return the secret.
        return msg