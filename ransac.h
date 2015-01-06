#ifndef RANSAC_H
#define RANSAC_H

#include <vector>
#include <math.h>
#include <combinatorics_utils.h>

namespace ransac {

template<class _Data>
class Model {
protected:
    typedef typename std::vector<_Data> DataSet;

public:
    /*!
      \return a measurment of how well the model fits to data
      */
    virtual double fit_to_vector(const std::vector<_Data> & data) = 0;

    /*!
      \return the distance between the model and a data point
    */
    virtual double distance_to_current(const _Data & query) const = 0;
}; // end class Model

////////////////////////////////////////////////////////////////////////////////

/*! input:
        \arg data - a set of observations
        \arg model - a model that can be fitted to data
        \arg n - the minimum number of data required to fit the model
        \arg k - the number of iterations performed by the algorithm
        \arg t - a threshold value for determining when a datum fits a model
        \arg d - the number of close data values required to assert that a model fits well to data
    output:
        \arg best_model - model parameters which best fit the data (or nil if no good model is found)
        \arg best_consensus_set - data point from which this model has been estimated
        \arg best_error - the error of this model relative to the data
    */
template<class _Data, class _Model>
static void ransac(const std::vector<_Data> & data,
                   const _Model & model,
                   const int minimium_data_for_model,
                   const int nb_iters,
                   const double thresh,
                   const size_t minimum_close_pts,
                   _Model & best_model,
                   std::vector<_Data> & best_consensus_set,
                   double & best_error
                   ) {
    //iterations := 0
    int iterations = 0;

    //best_model := nil

    //best_consensus_set := nil

    //best_error := infinity
    best_error = INFINITY;

    //while iterations < k
    while (iterations < nb_iters) {
       // printf("\niteration %i", iterations);

        //    maybe_inliers := n randomly selected values from data
        std::vector<_Data> maybe_inliers;
        // generate a random combination
        combinatorics_utils::UnorderedCombination maybe_inliers_combi;
        combinatorics_utils::combination_random(
                    maybe_inliers_combi, minimium_data_for_model, data.size());
        combinatorics_utils::select_random_sample_with_comb(
                    data, maybe_inliers, maybe_inliers_combi);

        //    maybe_model := model parameters fitted to maybe_inliers
        _Model maybe_model = model;
        maybe_model.fit_to_vector(maybe_inliers);

        //    consensus_set := maybe_inliers
        std::vector<_Data> consensus_set = maybe_inliers;

        //    for every point in data not in maybe_inliers
        for (unsigned int data_idx = 0; data_idx < data.size(); ++data_idx) {
            if (maybe_inliers_combi.count(data_idx) != 0)
                continue;

            const _Data* point = &data.at(data_idx);

            //    if point fits maybe_model with an error smaller than t
            if (maybe_model.distance_to_current(*point) < thresh) {
                //        add point to consensus_set
                consensus_set.push_back(*point);
            }
        }

        // if the number of elements in consensus_set is > d
        // (this implies that we may have found a good model,
        // now test how good it is)
        if (consensus_set.size() > minimum_close_pts) {
           // printf("\nit %i: foud a potential good model", iterations);

            // this_model := model parameters fitted to
            // all points in consensus_set
            // this_error := a measure of how well
            // this_model fits these points
            _Model this_model = model;
            double this_error = this_model.fit_to_vector(consensus_set);


            // if this_error < best_error
            // (we have found a model which is better than any of the previous ones,
            // keep it until a better one is found)
            if (this_error < best_error) {
               // printf("\nit %i: foud a really better model", iterations);

                // best_model := this_model
                best_model = this_model;

                // best_consensus_set := consensus_set
                best_consensus_set = consensus_set;

                // best_error := this_error
                best_error = this_error;
            } // end this_error < best_error
        } // end size > d

        // increment iterations
        ++iterations;
    } // end while iterations < k

    //return
} // end ransac()

} // end namespace ransac

#endif // RANSAC_H
