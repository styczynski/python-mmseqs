//
// Created by mad on 2019-02-27.
//

#ifndef MMSEQS_KMERMARKOVSCORE_H
#define MMSEQS_KMERMARKOVSCORE_H

#include <Indexer.h>

namespace MarkovScores {
static const int MARKOV_ORDER = 4;
static const float MEDIAN_SCORE = 2.025;

// trained matrix for std_plus_human
float markov5Scores[] = {
    1.484, 2.464, 1.910, 2.357, 1.753, 1.989, 2.101, 2.196, 1.843, 2.137, 1.849,
    2.210, 1.715, 1.858, 2.239, 2.269, 1.714, 2.389, 2.082, 1.899, 1.761, 2.257,
    2.131, 1.902, 2.509, 2.157, 1.819, 1.659, 2.117, 1.679, 2.318, 1.963, 1.591,
    2.316, 1.697, 2.653, 1.670, 2.094, 2.153, 2.142, 2.074, 2.257, 1.560, 2.223,
    1.823, 1.975, 2.281, 1.959, 1.454, 2.493, 1.942, 2.343, 1.997, 1.951, 2.092,
    1.963, 2.052, 2.207, 1.771, 2.004, 1.977, 1.681, 2.038, 2.392, 1.651, 2.120,
    2.078, 2.219, 2.158, 1.451, 2.394, 2.186, 2.369, 1.672, 1.908, 2.144, 2.163,
    1.533, 2.524, 1.958, 2.155, 2.033, 2.101, 1.747, 1.962, 2.221, 2.284, 1.628,
    2.925, 2.116, 1.871, 1.457, 2.233, 1.641, 2.393, 1.855, 1.833, 1.857, 1.676,
    2.941, 1.873, 1.963, 2.144, 2.034, 2.552, 1.751, 1.626, 2.262, 2.131, 1.823,
    2.456, 1.705, 2.036, 1.895, 1.817, 2.296, 2.499, 1.544, 2.226, 1.907, 2.549,
    1.722, 1.867, 1.989, 2.589, 1.288, 2.073, 2.421, 1.391, 2.510, 1.810, 2.662,
    1.880, 1.828, 2.088, 2.242, 2.149, 2.045, 1.545, 2.398, 1.824, 1.798, 2.093,
    2.353, 1.816, 2.228, 1.978, 2.007, 1.815, 2.239, 2.108, 1.878, 2.652, 2.079,
    1.697, 1.758, 1.999, 1.664, 2.332, 2.085, 1.643, 2.252, 1.632, 2.764, 1.766,
    2.023, 1.958, 2.304, 2.192, 2.137, 1.445, 2.424, 1.804, 1.830, 2.208, 2.212,
    1.635, 2.247, 1.765, 2.530, 2.226, 1.727, 2.082, 2.011, 2.256, 1.972, 1.669,
    2.176, 2.096, 1.574, 2.024, 2.438, 1.554, 2.395, 2.272, 1.931, 1.996, 1.795,
    2.453, 1.845, 2.099, 1.956, 2.071, 1.884, 1.823, 1.784, 2.376, 2.093, 2.010,
    2.209, 2.278, 1.602, 1.921, 2.217, 2.381, 1.605, 2.743, 2.106, 1.887, 1.524,
    2.368, 1.470, 2.425, 1.948, 1.717, 2.170, 1.917, 2.261, 1.893, 2.097, 2.266,
    1.790, 2.450, 1.908, 1.724, 2.013, 1.872, 1.916, 2.433, 1.854, 1.630, 2.400,
    1.979, 2.097, 2.294, 1.835, 2.404, 1.614, 2.339, 1.955, 2.054, 1.719, 2.195,
    1.495, 2.190, 2.267, 1.527, 2.335, 2.085, 2.191, 1.775, 2.079, 2.256, 1.934,
    1.964, 2.108, 1.931, 2.003, 1.994, 1.783, 2.456, 1.856, 1.868, 2.212, 1.958,
    1.984, 1.745, 2.458, 2.355, 1.625, 2.532, 2.275, 1.673, 1.705, 2.155, 1.567,
    2.487, 1.946, 1.718, 2.319, 1.760, 2.319, 1.738, 2.291, 2.378, 1.720, 2.174,
    2.237, 1.641, 2.026, 1.941, 1.880, 2.411, 1.838, 1.666, 2.270, 1.748, 2.474,
    1.771, 2.177, 2.414, 1.745, 2.208, 2.301, 1.579, 2.025, 2.068, 1.605, 2.092,
    2.335, 1.823, 2.020, 2.045, 2.130, 2.337, 1.471, 2.518, 1.906, 2.578, 1.654,
    1.967, 1.949, 2.393, 1.348, 2.561, 2.016, 2.286, 2.145, 2.127, 1.557, 1.883,
    2.338, 2.286, 1.616, 3.076, 2.091, 1.920, 1.386, 2.336, 1.556, 2.663, 1.718,
    1.926, 1.791, 1.668, 2.909, 2.214, 2.033, 2.415, 1.504, 2.907, 1.692, 1.846,
    1.841, 2.201, 1.575, 2.572, 1.843, 2.153, 1.834, 1.955, 2.079, 2.433, 1.557,
    2.463, 1.768, 2.637, 1.609, 2.002, 1.933, 2.431, 1.328, 2.292, 2.237, 1.360,
    2.495, 1.878, 2.635, 1.896, 1.934, 2.117, 2.065, 2.196, 1.983, 1.709, 2.166,
    1.775, 1.914, 2.065, 2.298, 1.742, 2.218, 1.955, 2.132, 1.759, 2.272, 1.937,
    2.081, 2.436, 2.189, 1.558, 1.963, 1.949, 1.775, 2.489, 1.886, 1.743, 2.183,
    1.627, 2.668, 1.846, 2.131, 2.029, 2.009, 2.148, 2.181, 1.567, 2.209, 1.887,
    1.784, 2.211, 2.163, 1.766, 2.109, 1.684, 2.616, 2.101, 2.044, 1.877, 1.988,
    2.279, 2.057, 1.530, 2.269, 2.059, 1.571, 2.146, 2.340, 2.014, 1.916, 2.426,
    1.731, 2.190, 1.869, 2.948, 1.406, 2.495, 1.675, 2.444, 1.618, 2.303, 1.434,
    2.865, 1.786, 2.143, 2.122, 2.144, 1.655, 2.070, 2.357, 2.619, 1.308, 3.134,
    2.073, 2.142, 1.245, 2.288, 1.400, 2.749, 1.902, 1.956, 1.997, 1.958, 2.093,
    2.067, 2.305, 2.708, 1.301, 2.583, 1.869, 1.906, 1.774, 2.086, 1.626, 2.730,
    1.787, 1.857, 2.127, 1.788, 2.285, 2.159, 2.028, 2.729, 1.396, 2.490, 1.960,
    2.300, 1.467, 2.368, 1.355, 2.231, 2.304, 1.437, 2.665, 1.951, 2.222, 1.836,
    2.054, 2.042, 2.082, 1.878, 2.304, 1.758, 2.121, 1.934, 1.886, 2.021, 2.175,
    1.704, 2.454, 1.983, 1.957, 1.852, 2.388, 2.086, 1.755, 2.357, 2.335, 1.637,
    1.811, 2.081, 1.721, 2.329, 1.936, 1.625, 2.546, 1.631, 2.460, 1.811, 2.214,
    1.995, 2.008, 2.139, 2.346, 1.484, 2.193, 1.837, 2.042, 2.149, 1.989, 1.523,
    2.610, 1.888, 2.197, 2.048, 1.948, 2.006, 2.000, 1.991, 2.292, 1.695, 2.087,
    2.023, 1.844, 1.866, 2.314, 1.804, 2.087, 1.960, 2.176, 2.496, 1.400, 2.335,
    2.026, 2.516, 1.577, 1.928, 2.138, 2.446, 1.429, 2.509, 1.893, 2.384, 2.168,
    2.055, 1.535, 1.955, 2.084, 2.271, 1.742, 3.152, 1.939, 1.884, 1.491, 2.582,
    1.576, 2.540, 1.619, 1.858, 1.880, 1.646, 2.912, 2.216, 1.829, 2.061, 1.924,
    2.939, 1.633, 1.564, 2.259, 2.310, 1.695, 2.351, 1.769, 2.252, 1.773, 1.881,
    2.145, 2.887, 1.394, 2.456, 1.727, 2.692, 1.454, 2.063, 2.054, 2.933, 1.106,
    2.217, 2.400, 1.433, 2.474, 1.754, 2.707, 1.931, 1.826, 1.967, 2.322, 2.221,
    2.012, 1.510, 2.424, 1.831, 1.800, 2.041, 2.406, 1.852, 2.318, 1.822, 2.062,
    1.749, 2.303, 1.912, 2.094, 2.568, 2.183, 1.491, 1.969, 2.079, 1.732, 2.249,
    1.989, 1.692, 2.358, 1.589, 2.617, 1.761, 2.088, 1.804, 2.447, 2.151, 2.187,
    1.486, 2.336, 1.887, 1.892, 2.014, 2.234, 1.702, 2.489, 1.644, 2.361, 2.188,
    1.736, 1.963, 2.160, 2.202, 2.097, 1.526, 2.309, 2.222, 1.604, 1.900, 2.406,
    1.665, 2.290, 2.234, 1.901, 2.285, 1.666, 2.571, 1.683, 2.365, 1.967, 2.049,
    1.697, 2.194, 1.534, 2.319, 2.085, 2.008, 2.269, 2.207, 1.611, 2.153, 2.206,
    2.585, 1.352, 3.129, 2.086, 2.151, 1.234, 2.707, 1.324, 2.599, 1.825, 1.820,
    2.204, 1.846, 2.175, 2.134, 2.104, 2.340, 1.548, 2.601, 1.928, 1.731, 1.883,
    2.153, 1.833, 2.280, 1.793, 1.685, 2.347, 1.909, 2.144, 2.322, 1.924, 2.403,
    1.525, 2.521, 2.037, 2.098, 1.521, 2.207, 1.519, 2.204, 2.202, 1.527, 2.375,
    1.948, 2.314, 1.752, 2.226, 2.093, 1.972, 1.961, 2.199, 1.800, 2.070, 1.755,
    1.878, 2.340, 2.094, 1.821, 2.348, 1.923, 1.961, 1.813, 2.497, 2.145, 1.679,
    2.492, 2.294, 1.639, 1.750, 1.852, 1.812, 2.519, 1.923, 1.636, 2.343, 1.740,
    2.461, 1.859, 2.294, 2.311, 1.649, 2.280, 2.316, 1.497, 2.065, 1.816, 1.873,
    2.376, 1.998, 1.592, 2.426, 1.765, 2.412, 1.832, 2.223, 2.211, 1.791, 2.192,
    2.312, 1.547, 2.073, 1.845, 1.750, 2.103, 2.384, 1.833, 2.049, 2.023, 2.110,
    2.386, 1.535, 2.572, 1.759, 2.652, 1.551, 1.989, 2.013, 2.350, 1.492, 2.690,
    1.769, 2.274, 2.181, 2.057, 1.589, 2.000, 2.406, 2.318, 1.471, 3.267, 2.138,
    1.877, 1.334, 2.114, 1.595, 2.594, 1.877, 1.847, 1.824, 1.648, 3.054, 1.984,
    2.281, 2.440, 1.485, 2.854, 1.732, 1.674, 2.015, 2.194, 1.758, 2.607, 1.637,
    2.129, 1.945, 1.767, 2.198, 2.311, 1.819, 2.361, 1.642, 2.667, 1.745, 1.961,
    1.799, 2.589, 1.248, 2.277, 2.276, 1.510, 2.400, 1.782, 2.568, 1.960, 1.977,
    2.064, 2.000, 2.397, 1.982, 1.499, 2.299, 1.877, 1.789, 2.092, 2.295, 2.092,
    2.174, 1.893, 1.866, 1.812, 2.416, 2.090, 1.771, 2.797, 2.142, 1.652, 1.683,
    1.882, 1.832, 2.399, 1.953, 1.682, 2.197, 1.657, 2.708, 1.884, 2.163, 2.138,
    1.844, 2.404, 2.083, 1.498, 2.177, 1.882, 1.809, 2.204, 2.144, 1.796, 2.196,
    1.675, 2.470, 2.227, 1.925, 1.947, 1.923, 2.373, 2.034, 1.622, 2.071, 2.153,
    1.643, 1.908, 2.408, 1.696, 2.285, 2.196, 1.900, 1.919, 2.201, 2.427, 1.591,
    2.301, 1.943, 2.006, 1.796, 1.834, 1.863, 2.429, 1.949, 1.875, 2.279, 2.016,
    1.867, 1.951, 2.556, 2.325, 1.428, 2.814, 2.247, 1.875, 1.417, 2.084, 1.708,
    2.516, 1.819, 1.743, 2.258, 1.832, 2.242, 1.772, 2.644, 2.386, 1.490, 2.503,
    2.089, 1.651, 1.888, 1.927, 1.870, 2.470, 1.822, 1.726, 2.432, 1.811, 2.137,
    2.066, 2.207, 2.395, 1.496, 2.436, 2.100, 1.971, 1.613, 2.020, 1.601, 2.142,
    2.341};
}  // namespace MarkovScores
class MarkovKmerScore {
 public:
  static float scoreKmer(const unsigned char* kmer, unsigned char kmerSize) {
    float totalSocore = 0.0;
    for (int pos = 0; pos < kmerSize - MarkovScores::MARKOV_ORDER; pos++) {
      size_t lookupIdx =
          Indexer::computeKmerIdx(&kmer[pos], MarkovScores::MARKOV_ORDER + 1);
      totalSocore += MarkovScores::markov5Scores[lookupIdx];
    }
    return totalSocore;
  }

  static int adjustedLength(const unsigned char* kmer, unsigned char kmerSize,
                            float minScoreThr) {
    float totalSocore = 0.0;
    int pos = 0;
    while (totalSocore < minScoreThr &&
           pos < kmerSize - MarkovScores::MARKOV_ORDER) {
      size_t lookupIdx =
          Indexer::computeKmerIdx(&kmer[pos], MarkovScores::MARKOV_ORDER + 1);
      float score = MarkovScores::markov5Scores[lookupIdx];
      totalSocore += score;
      pos++;
    }
    return pos + MarkovScores::MARKOV_ORDER;
  }
};

#endif  // MMSEQS_KMERMARKOVSCORE_H
