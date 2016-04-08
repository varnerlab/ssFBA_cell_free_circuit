# Script for the analysis of stocichiometric analysis -

function find_source_nodes(stocichiometric_matrix::Array{Float64,2})

  # Find nodes that are act as sources in the network -
  (number_of_rows,number_of_cols) = size(stocichiometric_matrix);

  # return array -
  source_node_array = zeros(number_of_rows);

  # Go thru each row, is there are only negative entries?
  for metabolite_index in collect(1:number_of_rows)

    # Get the row -
    stm_row = stocichiometric_matrix[metabolite_index,:];

    # How many netaive elements are there?
    idx_negative = find(stm_row.<0);

    # How many postive elements are there?
    idx_positive = find(stm_row.>0);

    # check -
    if (length(idx_positive) == 0 && length(idx_negative) != 0)
      source_node_array[metabolite_index] = 1;
    end
  end

  return source_node_array;
end

function find_sink_nodes(stocichiometric_matrix::Array{Float64,2})

  # Find nodes that are act as sources in the network -
  (number_of_rows,number_of_cols) = size(stocichiometric_matrix);

  # return array -
  sink_node_array = zeros(number_of_rows);

  # Go thru each row, is there are only negative entries?
  for metabolite_index in collect(1:number_of_rows)

    # Get the row -
    stm_row = stocichiometric_matrix[metabolite_index,:];

    # How many netaive elements are there?
    idx_negative = find(stm_row.<0);

    # How many postive elements are there?
    idx_positive = find(stm_row.>0);

    # check -
    if (length(idx_positive) !=0 && length(idx_negative) == 0)
      sink_node_array[metabolite_index] = 1;
    end
  end

  return sink_node_array;
end
