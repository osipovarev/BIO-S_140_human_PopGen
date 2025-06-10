### The following code simulates genetic drift in populations of different sizes

# Function to simulate genetic drift
simulate_genetic_drift <- function(population_size, generations, allele_freq) {
  # Initialize population with initial allele frequencies
  population <- rep(0, population_size)
  population[1:round(allele_freq * population_size)] <- 1
  
  # Simulate generations
  for (gen in 1:generations) {
    # Calculate allele frequencies in the current generation
    allele_freq <- sum(population) / population_size
    
    # Plot the data point
    points(gen, allele_freq, type="o", pch=20, col="red")
    points(gen, 1 - allele_freq, type="o", pch=20, col='blue')
    
    # Sample the next generation based on the current allele frequencies
    offspring <- sample(population, size = population_size, replace = TRUE)
    
    # Update population with the offspring
    population <- offspring
  }
  
  # Return final allele frequency
  final_allele_freq <- sum(population) / population_size
  # return(final_allele_freq)
}


# Set the parameters for the simulation!!
population_size <- 10
generations <- 50
allele_freq <- 0.5


## Initiate plotting
plot(0, allele_freq, 
     xlim = c(1, generations),
     ylim = c(0, 1),
     main = "genetic drift",
     xlab = 'number of generations',
     ylab = 'allele frquency', 
     frame.plot = FALSE)


# Run the simulation!!
simulate_genetic_drift(population_size, generations, allele_freq)



#simulate_genetic_drift(population_size, generations, allele_freq)

# Print the final allele frequency
# print(paste("Final allele frequency:", final_freq))
  
