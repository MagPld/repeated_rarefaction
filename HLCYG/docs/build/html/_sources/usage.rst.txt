Usage
=====

There are two main functions to the package. The first one repeated_rarefaction (instert reference link here) performs repeated rarefaction, caluclates an ordination plot with NMDS and plots the resulting ordination. The second allows for testing a range of thresholds with repeated rarefaction, allowing the user to see the effect of lowering the rarefaction threshold as well as the use of different amount of repeats.

.. function:: repeated_rarefaction(physeq, repeats, threshold, method, sample_id, colorb, shapeb, cloud = FALSE, ellipse = TRUE)
   
   Perform one iteration of repeated rarefaction and produce ordination plot

   :param physeq: Object containing the count data and information about it.
   :type physeq: phyloseq object
   :param repeats: Indicates the amount of repeats run. If repeats = 1 only one run and therefore performed without repeats.
   :type repeats: integer
   :param threshold: The threshold value for the rarefaction.
   :type threshold: integer
   :param method: Method used for the ordination caluclation. Currently only "NMDS" supported.
   :type method: string
   :param sample_id: Name of the column with sample IDs.
   :type sample_id: string
   :param colorb: Name of the column with data to colour the graph by.
   :type colorb: string
   :param shapeb: Name of the column with data to shape the points on the graph by.
   :type shapeb: string
   :param cloud: Aesthetic setting for the graph. Default is FALSE, showing only the median point for each sample. TRUE shows datapoints for all repeats.
   :type cloud: boolean.
   :param ellipse: Aesthetic setting for the graph. Default is TRUE, drawing the 95% confidence ellipses around each sample.
   :type ellipse: boolean.

.. function:: test_threshold(physeq, repeats, t_min, t_max, t_step, method = "NMDS", sample_id, groupb)

   Test a range of threshold values for repeated rarefaction.

   :param physeq: Object containing the count data and information about it.
   :type physeq: phyloseq object
   :param repeats: Indicates the amount of repeats run. If repeats = 1 only one run and therefore performed without repeats.
   :type repeats: integer
   :param t_min: Minimum value for the threshold testing range.
   :type t_min: integer
   :param t_max: Maximum value for the threshold testing range.
   :type t_max: integer
   :param t_step: Step value for the threshold testing range. A value between 0 and 1 will cause the same threshold to be tested more than once.
   :type t_step: float
   :param method: Method used for the ordination caluclation. Currently only "NMDS" supported.
   :type method: string
   :param sample_id: Name of the column with sample IDs.
   :type sample_id: string
   :param groupb: Name of the column with data to group the points by.
   :type groupb: string