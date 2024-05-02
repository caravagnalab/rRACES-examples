# adding mutations
# add mutations

# Mutation generation ####
library(rRACES)
# m_engine <- build_mutation_engine(setup_code = "GRCh38", context_sampling = 50)

# m_engine <- build_mutation_engine(reference_src = "./GRCh38/reference.fasta", SBS_src = "./GRCh38/SBS.txt", 
#                                   drivers_src = "./GRCh38/drivers.txt", passenger_CNAs_src = "./GRCh38/passenger_CNAs.txt", 
#                                   germline_src = "./GRCh38/germline_data/", directory = "./GRCh38/")

m_engine <- build_mutation_engine(setup_code = "GRCh38")

# chr7	140753332	140753332	T	A	SNV	CLLE	BRAF	BRAF V600E
m_engine$add_mutant("Clone 1", c(SNV = 5e-9), c(SNV("7", 140753332, "A"))) 


m_engine$add_mutant(mutant_name = "Clone 2", 
                    passenger_rates = c(SNV = 3e-8, CNA = 1e-11),
                    driver_SNVs = SNV("3",  179218303, "A"))

m_engine$add_mutant(mutant_name = "Clone 3", 
                    passenger_rates = c(SNV = 3e-5, CNA = 1e-11), # modelled as an hypermutant, so with higher passenger_rates
                    driver_SNVs = SNV("2", 47799065, "A"))

m_engine$add_exposure(c(SBS1 = 0.5,SBS5 = 0.5))
m_engine$add_exposure(22, c(SBS1 = 0.3,SBS5 = 0.3, SBS6 = 0.4))

