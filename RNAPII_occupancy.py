import numpy as np
import matplotlib.pyplot as plt


#params: list of parameters [K_4, K_2, K_-2, K_-3]
#n_genes: number of genes
#K_1: recruitment rate of RNAPII
#K_1r: dissociation rate of RNAPII at UAS
#K_5: initiation rate
#K_6: elongation rate
#K_7: termination rate
def STM(params, n_genes=1000, K_1 = 0.002, K_1r =0.003, K_5=0.1, K_6=0.1388, K_7=0.0325):
    
    gene_history = [] #container for gene states
    gene_times = [] #container for times
    
    
    #  [recruitment of GTFs, recruitment of RNAPII, dissociation at UAS, UAS->Pro, Pro->UAS, dissociation at promoter, Initiation, elongation (x10), termination]
    transition_rates = [params[0], K_1, K_1r, params[2], params[3], params[1], K_5] + [K_6]*10 + [K_7] 
   
    T = 1000 #duration
    
    for gene in range(n_genes):
        
        states = [0]*(len(transition_rates)-4)
        t = 0
        times = [t]
        states_history = []
    
        while t < T: #solve using Gillespie alg
            states_history.append(np.copy(states))
            R0 = transition_rates[0]*(states[0]==0)*(states[2]==1) #late GTF comes in only if there is a spot available
            R1 = R0 + transition_rates[1]*(states[1]==0) #polii come in if theres a spot available
            R2 = R1 + transition_rates[2]*(states[1]==1) #polii go out if its already there
            R3 = R2 + transition_rates[3]*(states[1]==1)*(states[2]==0) #polii go to promoter if it was at UAS
            R4 = R3 + transition_rates[4]*(states[1]==0)*(states[2]==1) #polii go back to enhancer if theres not one there
            R5 = R4 + transition_rates[5]*(states[2]==1) #polii go out if its at promoter
            R6 = R5 + transition_rates[6]*(states[2]==1)*(states[0]==1) #starting elongation if there are late GTFs
            Rs = [R0, R1, R2, R3, R4, R5, R6]
            for i in np.arange(len(Rs),len(transition_rates)-1): #elongating
                Rs.append(Rs[-1] + transition_rates[i]*(states[i-3]==0)*(states[i-4]==1))
            Rs.append(Rs[-1] + transition_rates[-1]*(states[-1]==1)) #stuck at the terminator for some time
            Rs = np.array(Rs)
            R = Rs[-1]
    
            # Calculate time to next reaction
            u1 = np.random.random()
            tau = 1/R * np.log(1/u1)
    
            # Select which reaction to occur
            u2 = np.random.random()
            j = np.argwhere(Rs > u2*R)[0][0]
    
            if j == 0: #GTF comes in, allowing initiation
                states[0] = 1
            elif j == 1: #pol2 comes in 
                states[1] = 1
            elif j == 2: #pol2 goes out
                states[1] = 0
            elif j == 3: #movement to promoter
                states[2] = 1
                states[1] = 0
            elif j == 4: #movement back to enhancer
                states[1] = 1
                states[2] = 0
            elif j == 5: #pol2 goes out
                states[2] = 0
                states[0] = 0
            elif j == 6: #elongating starts
                states[3] = 1
                states[2] = 0
                states[0] = 0
            elif j == (len(Rs)-1): #leaves terminator
                states[-1] = 0
            else: #elongating
                states[j-3] = 1 #move one spot down
                states[j-4] = 0
    
            t += tau
            times.append(t)
            
        times[-1] = T
        states_history = np.array(states_history)
        
        
        #consider only final minute as in the experiments
        measurement_time = 60
        t_idx = np.argwhere(np.array(times)>=(T-measurement_time))[0][0]
        gene_history.append(states_history[t_idx-1:,:])
        gene_times.append(times[t_idx-1:])
        gene_times[gene][0] = T-measurement_time

    counts = 0.0*np.array(states)
    for gene in range(n_genes):
        dt = np.diff(gene_times[gene])
        counts += gene_history[gene].T@dt
    enh = counts[1]/n_genes
    pro = counts[2]/n_genes
    cds = sum(counts[3:13])/n_genes/10 #mean
    term= counts[-1]/n_genes
    data = np.array([enh,pro,cds,term])
    return data

#params: list of parameters [K_4, K_-3]
#n_genes: number of genes
#K_3: recruitment rate of RNAPII
#K_5: initiation rate
#K_6: elongation rate
#K_7: termination rate
def TFO(params, n_genes=1000, K_3 = 0.002, K_5=0.1, K_6 = 0.1388, K_7=0.037):
    
    gene_history = [] #container for gene states
    gene_times = [] #container for times

    #  [recruitment of GTFs, recruitment of RNAPII, dissociation at promoter, Initiation, elongation (x10), termination]
    transition_rates = [params[0], K_3, params[1], K_5] + [K_6]*10 + [K_7] 

    
    T = 1000 #duration
    
    for gene in range(n_genes):
        
        states = [0]*(len(transition_rates)-2)
        t = 0
        times = [t]
        states_history = []
    
        while t < T: #solve using Gillespie alg
            states_history.append(np.copy(states))
            
            R0 = transition_rates[0]*(states[0]==0)*(states[1]==1) #GTFs in only if its not already there
            R1 = R0 + transition_rates[1]*(states[1]==0) #polii come in only if there is a spot available 
            R2 = R1 + transition_rates[2]*(states[1]==1) #polii go out only if its already there
            R3 = R2 + transition_rates[3]*(states[1]==1)*(states[0]==1) #polii initiate only if polii and GTFs are there
        
            Rs = [R0, R1, R2, R3]
            
            n_start = len(Rs)
            for i in np.arange(n_start,len(transition_rates)-1): #elongating
                Rs.append(Rs[-1] + transition_rates[i]*(states[i-1]==0)*(states[i-2]==1))
                
            Rs.append(Rs[-1] + transition_rates[-1]*(states[-1]==1)) #stuck at the terminator for some time
            Rs = np.array(Rs)
            R = Rs[-1]
    
            # Calculate time to next reaction
            u1 = np.random.random()
            tau = 1/R * np.log(1/u1)
    
            # Select which reaction to occur
            u2 = np.random.random()
            j = np.argwhere(Rs > u2*R)[0][0]
    
            if j == 0: #TFIIH comes in
                states[0] = 1
            elif j == 1: #polii comes in
                states[1] = 1
            elif j == 2: #polii goes out
                states[1] = 0
            elif j == 3: #polii initiates
                states[0] = 0
                states[1] = 0
                states[2] = 1
            elif j == (len(Rs)-1): #leaves terminator
                states[-1] = 0
            else: #elongating
                states[j-1] = 1 #move one spot down
                states[j-2] = 0
    
            t += tau
            times.append(t)
            
        times[-1] = T
        states_history = np.array(states_history)
        
        #consider only final minute
        measurement_time = 60
        t_idx = np.argwhere(np.array(times)>=(T-measurement_time))[0][0]
        gene_history.append(states_history[t_idx-1:,:])
        gene_times.append(times[t_idx-1:])
        gene_times[gene][0] = T-measurement_time

    counts = 0.0*np.array(states)
    for gene in range(n_genes):
        dt = np.diff(gene_times[gene])
        counts += gene_history[gene].T@dt
    pro = counts[1]/n_genes
    cds = sum(counts[2:12])/n_genes/10 #mean
    term= counts[-1]/n_genes
    data = np.array([pro,cds,term])
    return data


def main():
    STM_example = STM([0.013, 0.025, 0.052, 0.001])
    TFO_example = TFO([0.022, 0.019])

    fig, axs = plt.subplots(1,2, figsize=(12,3))
    axs[0].bar(np.array([0,1,2,3]), STM_example)

    axs[0].set_xticks([0,1,2,3], labels=["UAS", "Pro", "CDS", "Term"])
    axs[1].bar([0,1,2], TFO_example)
    axs[1].set_xticks([0,1,2], labels=["Pro", "CDS", "Term"])
    axs[0].set_ylabel("Dwell time")
    axs[0].set_title("STM example")
    axs[1].set_title("TFO example")
    axs[0].set_ylim(0,5)
    axs[1].set_ylim(0,5)
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()