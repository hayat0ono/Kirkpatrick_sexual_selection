# 基本個体クラスの定義

class Individual:
    def __init__(self, t_gene=0, p_gene=0):
        """
        Args:
            t_gene (int): T gene (phenotype) - 0 or 1
            p_gene (int): P gene (preference) - 0 or 1
        """
        self.t_gene = t_gene
        self.p_gene = p_gene
        self.death_prob = 0.0
        self.year = 0
    
    def get_t(self):
        return self.t_gene
    
    def get_p(self):
        return self.p_gene
    
    def get_age(self):
        return self.year
    
    def set_death_prob(self, prob):
        self.death_prob = prob
    
    def age(self):
        self.year += 1 