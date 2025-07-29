# モデルクラスの定義

import random
import numpy as np
from male import Male
from female import Female

class Model:    
    def __init__(self, N=1000, male_female_ratio=0.5, t_gene_ratio=0.5, p_gene_ratio=0.5,
                 male_death_prob=0.0, s_ratio=0.2, female_death_prob=0.0,
                 a0_coeff=2.0, a1_coeff=3.0, end_time=100, lifetime=1, num_child=2,
                 mutation_rate=0.0):
        """
        Args:
            N (int): Total number of individuals
            male_female_ratio (float): Ratio of males to total population
            t_gene_ratio (float): Initial frequency of T1 gene in total population
            p_gene_ratio (float): Initial frequency of P1 gene in total population
            male_death_prob (float): Base death probability for males before mating
            s_ratio (float): Difference in death probability between T gene males (T1 males have (1-s_ratio) times survival rate compared to T0 males)
            female_death_prob (float): Death probability for females before mating
            a0_coeff (float): Selection coefficient in Kirkpatrick's theory for female preference
            a1_coeff (float): Selection coefficient in Kirkpatrick's theory for female preference
            end_time (int): Number of generations to run
            lifetime (int): Lifespan of individuals
            num_child (int): Number of offspring per mating
            mutation_rate (float): Mutation rate (0.0-1.0 range)
        """
        self.N = N
        self.male_female_ratio = male_female_ratio
        self.t_gene_ratio = t_gene_ratio
        self.p_gene_ratio = p_gene_ratio
        self.a0_coeff = a0_coeff
        self.a1_coeff = a1_coeff
        self.male_death_prob = male_death_prob
        self.female_death_prob = female_death_prob
        self.s_ratio = s_ratio
        self.end_time = end_time
        self.lifetime = lifetime
        self.num_child = num_child
        self.mutation_rate = mutation_rate
        
        # Internal variables
        self.t_male_death_prob = [0.0, 0.0] 
        self.breed_prob = [0.0, 0.0] 
        self.figures = np.zeros((2, 2, 2), dtype=int) 
        self.time = 0 
        
        # Individual lists
        self.t_male_list = [[], []]  
        self.female_list = [] 
        
        # Coefficients for gene frequency relationship at equilibrium in Kirkpatrick's theory
        self.v1 = 0.0
        self.v2 = 0.0
        
        # Initialize
        self._init_parameters()
    
    def _init_parameters(self):
        self.t_male_death_prob[0] = self.male_death_prob
        self.t_male_death_prob[1] = 1.0 - ((1.0 - self.male_death_prob) * (1.0 - self.s_ratio))
        
        self.v1 = (self.a0_coeff + self.s_ratio - 1) / (self.a0_coeff * self.a1_coeff - 1) / (1.0 - self.s_ratio)
        self.v2 = self.a1_coeff * (self.a0_coeff + self.s_ratio - 1) / (self.a0_coeff * self.a1_coeff - 1)
    
    def _apply_mutation(self, t_gene, p_gene):
        new_t_gene = t_gene
        new_p_gene = p_gene
        
        # T gene mutation
        if random.random() < self.mutation_rate:
            new_t_gene = 1 - t_gene
        
        # P gene mutation
        if random.random() < self.mutation_rate:
            new_p_gene = 1 - p_gene
        
        return new_t_gene, new_p_gene
    
    def _clear_lists(self):
        for i in range(2):
            self.t_male_list[i].clear()
        self.female_list.clear()
    
    def build_objects(self):
        self._clear_lists()
        self._init_parameters()
        self.time = 0
        
        # Initialize statistics
        self.figures.fill(0)
        
        # Create individuals
        for _ in range(self.N):
            # Determine T gene
            if random.random() < self.t_gene_ratio:
                t_gene = 1
            else:
                t_gene = 0
            
            # Determine P gene
            if random.random() < self.p_gene_ratio:
                p_gene = 1
            else:
                p_gene = 0
            
            # Apply mutation
            t_gene, p_gene = self._apply_mutation(t_gene, p_gene)
            
            # Determine sex
            if random.random() < self.male_female_ratio:
                # Create male individual
                a_male = Male(t_gene, p_gene)
                a_male.set_death_prob(self.t_male_death_prob[t_gene])
                self.t_male_list[t_gene].append(a_male)
                self._add_figures(0, t_gene, p_gene)
            else:
                # Create female individual
                a_female = Female(t_gene, p_gene)
                a_female.set_death_prob(self.female_death_prob)
                self.female_list.append(a_female)
                self._add_figures(1, t_gene, p_gene)
        
        print(f"Death Probability\nT0Male: {self.t_male_death_prob[0]}\tT1Male: {self.t_male_death_prob[1]}\tFemale: {self.female_death_prob}")
        print(f"Mutation Rate: {self.mutation_rate}")
        self._calc_breed_prob()
        self.show_ratio()
        self.show_count()
        self.time += 1
    
    def init_objects(self):
        self.build_objects()
    
    def count_t_male(self, t_gene):
        return len(self.t_male_list[t_gene])
    
    def count_female(self):
        return len(self.female_list)
    
    def count_t0_gene(self):
        num_t = 0
        for x in range(2):
            for y in range(2):
                num_t += self.figures[x][0][y]
        return num_t
    
    def count_t1_gene(self):
        num_t = 0
        for x in range(2):
            for y in range(2):
                num_t += self.figures[x][1][y]
        return num_t
    
    def count_p0_gene(self):
        num_p = 0
        for x in range(2):
            for y in range(2):
                num_p += self.figures[x][y][0]
        return num_p
    
    def count_p1_gene(self):
        num_p = 0
        for x in range(2):
            for y in range(2):
                num_p += self.figures[x][y][1]
        return num_p
    
    def count_total(self):
        return self.count_t_male(0) + self.count_t_male(1) + self.count_female()
    
    def get_t1_ratio(self):
        total = self.count_total()
        if total == 0:
            return 0.0
        return self.count_t1_gene() / total
    
    def get_p1_ratio(self):
        total = self.count_total()
        if total == 0:
            return 0.0
        return self.count_p1_gene() / total
    
    def check_to_stop(self):
        if self.time >= self.end_time:
            print("Timeout")
            return False
        elif self.count_t_male(0) == 0 and self.count_t_male(1) == 0:
            print("Male dies out")
            return False
        elif self.count_female() == 0:
            print("Female dies out")
            return False
        
        self.time += 1
        return True
    
    def show_count(self):
        print(f"\t (T0, T1) = ({self.count_t0_gene()}, {self.count_t1_gene()})\t"
              f"(P0, P1) = ({self.count_p0_gene()}, {self.count_p1_gene()})")
    
    def show_ratio(self):
        print(f"[Generation:{self.time}] T1: {self.get_t1_ratio():.4f}\tP1: {self.get_p1_ratio():.4f}")
    
    def _calc_breed_prob(self):
        t0_male = self.count_t_male(0)
        t1_male = self.count_t_male(1)
        
        if t1_male + self.a0_coeff * t0_male > 0:
            self.breed_prob[0] = (self.a0_coeff * t0_male) / (t1_male + self.a0_coeff * t0_male)
        else:
            self.breed_prob[0] = 0.0
            
        if t0_male + self.a1_coeff * t1_male > 0:
            self.breed_prob[1] = (self.a1_coeff * t1_male) / (t0_male + self.a1_coeff * t1_male)
        else:
            self.breed_prob[1] = 0.0
    
    def _coupling(self):
        # Shuffle lists
        random.shuffle(self.t_male_list[0])
        random.shuffle(self.t_male_list[1])
        random.shuffle(self.female_list)
        
        female_count = self.count_female()
        
        for _ in range(female_count):
            # Take female from front of list
            a_female = self.female_list.pop(0)
            
            if random.random() < self.breed_prob[a_female.get_p()]:
                # P0 mates with T0, P1 mates with T1
                self._select(a_female.get_t(), a_female.get_p(), a_female.get_p())
            else:
                # P0 mates with T1, P1 mates with T0
                self._select(a_female.get_t(), a_female.get_p(), (a_female.get_p() + 1) % 2)
            
            # Return female to end of list
            self.female_list.append(a_female)
    
    def _select(self, female_t, female_p, male_t):
        male_t_count = self.count_t_male(male_t)
        
        # If no T gene male exists for mating, no mating occurs
        if male_t_count <= 0:
            return
        
        for _ in range(male_t_count):
            # Take male from list
            a_male = self.t_male_list[male_t].pop(0)
            
            # If not yet mated, perform mating
            if not a_male.married:
                a_male.married = True
                for _ in range(self.num_child):
                    self._breed(a_male.get_t(), a_male.get_p(), female_t, female_p)
                # Return male to end of list
                self.t_male_list[male_t].append(a_male)
                break
            # Return male to end of list
            self.t_male_list[male_t].append(a_male)
    
    def _breed(self, male_t, male_p, female_t, female_p):
        # Determine offspring T gene
        if male_t == female_t:
            t_gene = male_t
        else:
            if random.random() < 0.5:
                t_gene = male_t
            else:
                t_gene = female_t
        
        # Determine offspring P gene
        if male_p == female_p:
            p_gene = male_p
        else:
            if random.random() < 0.5:
                p_gene = male_p
            else:
                p_gene = female_p
        
        # Apply mutation
        t_gene, p_gene = self._apply_mutation(t_gene, p_gene)
        
        # Determine sex
        if random.random() < 0.5:
            # Create male
            a_male = Male(t_gene, p_gene)
            a_male.set_death_prob(self.t_male_death_prob[t_gene])
            self.t_male_list[t_gene].append(a_male)
            self._add_figures(0, t_gene, p_gene)
        else:
            # Create female
            a_female = Female(t_gene, p_gene)
            a_female.set_death_prob(self.female_death_prob)
            self.female_list.append(a_female)
            self._add_figures(1, t_gene, p_gene)
    
    def _add_figures(self, sex, t_gene, p_gene):
        self.figures[sex][t_gene][p_gene] += 1
    
    def _remove_figures(self, sex, t_gene, p_gene):
        self.figures[sex][t_gene][p_gene] -= 1
    
    def _die(self):
        # Male death processing
        for x in range(2):
            # Get number of Tx gene males
            count = self.count_t_male(x)
            
            for _ in range(count):
                # Take individual from front of list
                a_male = self.t_male_list[x].pop(0)
                
                if (a_male.get_age() >= self.lifetime or 
                    random.random() < self.t_male_death_prob[x]):
                    # Death
                    self._remove_figures(0, x, a_male.get_p())
                else:
                    # Return individual to end of list
                    self.t_male_list[x].append(a_male)
        
        # Female death processing
        count = self.count_female()
        
        for _ in range(count):
            # Take individual from front of list
            a_female = self.female_list.pop(0)
            
            if (a_female.get_age() >= self.lifetime or 
                random.random() < self.female_death_prob):
                # Death
                self._remove_figures(1, a_female.get_t(), a_female.get_p())
            else:
                # Return individual to end of list
                self.female_list.append(a_female)
    
    def _next_age(self):
        # Males
        for x in range(2):
            count = self.count_t_male(x)
            
            for _ in range(count):
                a_male = self.t_male_list[x].pop(0)
                a_male.age()
                self.t_male_list[x].append(a_male)
        
        # Females
        count = self.count_female()
        
        for x in range(count):
            a_female = self.female_list.pop(0)
            a_female.age()
            self.female_list.append(a_female)
    
    def get_equilibrium_value(self):
        p1_ratio = self.get_p1_ratio()
        tmp = 1.0 - self.s_ratio
        
        if p1_ratio <= self.v1:
            return 0.0
        elif p1_ratio < self.v2:
            return ((self.a0_coeff * self.a1_coeff - 1) * tmp / (self.a0_coeff - tmp) * p1_ratio - 1) / (self.a1_coeff * tmp - 1)
        else:
            return 1.0
    
    def step(self):
        # Calculate female preference probability by gene
        self._calc_breed_prob()
        
        # Increase age by 1 and reset male married status to false
        self._next_age()
        
        # Mating
        self._coupling()
        
        # Individuals die with probability or when reaching lifetime
        self._die()
        
        # Display current situation
        self.show_ratio()
        self.show_count() 