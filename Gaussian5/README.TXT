In this example we have developed different scenarios for model fitting,
and they are split in different file (see below). The aim is to explore
how AMIS with INLA estimates the model under a number of scenarios.

* 02A_AMIS.R: AMIS with 'standard' setting (5000 + 10 * 1000 samples) and vague initial parameters for the sampling distribution.
* 02B_AMIS.R: AMIS with 'standard' setting (5000 + 10 * 1000 samples) and parameters informed from the data for the sampling distribution.
* 02C_AMIS.R: AMIS with  1000 + 10 * 1000 samples and vague initial parameters for the sampling distribution.
* 02D_AMIS.R: AMIS with 1000 + 10 * 1000 samples and parameters informed from the data for the sampling distribution.
* 02B_AMIS-10.R:  AMIS with 'standard' setting (5000 + 10 * 1000 samples) and vague initial parameters for the sampling distribution but with a higher variance.
* 02B_AMIS-55K.R: AMIS with 5000 + 10 * 5000 samples and parameters informed from the data for the sampling distribution.

