# 雄個体クラスの定義

from individual import Individual

class Male(Individual):
    def __init__(self, t_gene=0, p_gene=0):
        """        
        Args:
            t_gene (int): T gene (phenotype) - 0 or 1
            p_gene (int): P gene (preference) - 0 or 1
        """
        super().__init__(t_gene, p_gene)
        self.married = False
    
    def age(self):
        super().age()
        self.married = False 