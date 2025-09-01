#include "mex.h"
#include <algorithm>  // for std::max, std::min

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Input check
    if (nrhs != 4)
        mexErrMsgIdAndTxt("clip_spikes:nrhs", "Four inputs required: dataRAW (single), spike_row, spike_col, half_width.");

    if (nlhs != 1)
        mexErrMsgIdAndTxt("clip_spikes:nlhs", "One output required: data_copy.");

    // Validate dataRAW
    if (!mxIsSingle(prhs[0]))
        mexErrMsgIdAndTxt("clip_spikes:type", "Input dataRAW must be of type single.");

    // Validate spike_row and spike_col
    if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
        mexErrMsgIdAndTxt("clip_spikes:type", "spike_row and spike_col must be of type double.");

    // Validate half_width
    if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1)
        mexErrMsgIdAndTxt("clip_spikes:type", "half_width must be a scalar of type double.");

    int half_width = static_cast<int>(mxGetScalar(prhs[3]));
    if (half_width < 0)
        mexErrMsgIdAndTxt("clip_spikes:value", "half_width must be non-negative.");

    // Input data
    float *dataRAW = static_cast<float *>(mxGetData(prhs[0]));
    mwSize num_rows = mxGetM(prhs[0]);
    mwSize num_cols = mxGetN(prhs[0]);

    double *spike_row = mxGetPr(prhs[1]);
    double *spike_col = mxGetPr(prhs[2]);

    mwSize num_spikes = std::min(mxGetNumberOfElements(prhs[1]), mxGetNumberOfElements(prhs[2]));

    // Output array (copy of dataRAW)
    plhs[0] = mxDuplicateArray(prhs[0]);
    float *data_copy = static_cast<float *>(mxGetData(plhs[0]));

    // Loop through spikes
    for (mwSize i = 0; i < num_spikes; ++i) {
        int r = static_cast<int>(spike_row[i]) - 1; // Convert to 0-based indexing
        int c = static_cast<int>(spike_col[i]) - 1;

        if (r < 0 || r >= static_cast<int>(num_rows) || c < 0 || c >= static_cast<int>(num_cols))
            continue;

        int r_start = std::max(r - half_width, 0);
        int r_end   = std::min(r + half_width, static_cast<int>(num_rows) - 1);

        float val_start = data_copy[r_start + c * num_rows];
        float val_end   = data_copy[r_end   + c * num_rows];
        int span = r_end - r_start;

        float slope = span > 0 ? (val_end - val_start) / span : 0.0f;

        for (int row = r_start; row <= r_end; ++row) {
            float interp_val = val_start + (row - r_start) * slope;
            data_copy[row + c * num_rows] = interp_val;
        }
    }
}
