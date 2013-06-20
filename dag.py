from pprint import pprint
import sys
from operator import itemgetter
import csv

COL = 0
ROW = 1

class Node(object):
    """docstring for Node"""
    def __init__(self, weight=0, parent=[], match=0):
        self.weight = weight
        # parent is (x,y) coord
        self.parent = parent
        self.match = match

class contig_node(Node):
    """this represents a contig in a special DAG"""
    def __init__(self, contig_path, weight, parent_node=[]):
        self.weight = weight
        self.parent_node = parent_node
        self.contig_path = contig_path

class DAG(object):
    """docstring for DAG"""
    def __init__(self, list1, list2):
        self.list1 = list1
        self.list2 = list2

        self.num_rows = len(self.list1)
        self.num_cols = len(self.list2)

        self.dag_array = None
        self.multi_paged_array = None
        self.list_of_contigs = None
        self.NON_MATCH_CHAR = '#'
        self.match_weight = 10
        self.miss_weight = 30
        self.indel = 20
        self.word_size = self.match_weight*5

    def compute_all_long_contigs(self):
        self.compare ()
        self.make_long_contigs_list()

    def make_long_contigs_list(self):
        list_of_contigs = [[]]
        for col in range(self.num_cols -1, 0, -1):
            for row in range(self.num_cols -1, 0, -1):
                self.get_contig_if_missing(col, row, list_of_contigs)
        self.list_of_contigs = [contigs for contigs in list_of_contigs if contigs]

    def remove_postfix_mismatches_from_contigs(self, contigs):
        if len(contigs) < 2:
            return
        cur_node = 0
        prev_node = self.get(contigs[0]).weight
        to_del = []
        for count, contig in enumerate(contigs[1:]):
            cur_node = self.get(contig).weight
            if prev_node < cur_node:
                to_del.append(count)
            else:
                break
            prev_node = cur_node
        for count in to_del:
            del contigs[0]

    def get_contig_if_missing(self, col, row, list_of_contigs):
        if self.dag_array[col][row].weight < self.word_size:
            return
        if self.check_if_exists(self.dag_array[col][row], list_of_contigs):
            return
        # it's new and big enough
        coord_path = self.get_path_from_coord((col, row))
        if coord_path:
          self.remove_postfix_mismatches_from_contigs(coord_path)
        for coord in coord_path:
            if self.check_if_exists(self.get(coord), list_of_contigs):
                return
        list_of_contigs.append(coord_path)
    
    def check_if_exists(self, the_node, list_of_nodes):
        for a_node_chain in list_of_nodes:
            if the_node.parent in [self.get(a_node).parent for a_node in a_node_chain]:
                return True
        return False

    def get_path_from_coord(self, node_coord):
        node_path_list = [] 
        cur_node = self.get(node_coord)
        coord = node_coord[:]
        while cur_node.weight > 0:
            node_path_list.append(coord)
            coord = cur_node.parent['node']
            cur_node = self.get(coord)
        return node_path_list

    def compare(self):
        # Create 2d array
        self.dag_array = [[None for col in range(self.num_cols)] for row in range(self.num_rows)]
        for col in range(0, self.num_cols):
            for row in range(0, self.num_rows):
                self.dag_array[col][row] = self.get_weighted_node(col, row)

    def get_weighted_node(self, col, row):
        match = self.list1[col] == self.list2[row]

        # First item, no parent
        if row == col == 0:
            if match:
                return Node(self.match_weight, None, match)
            else:
                return Node(0, None, match)

        # If we get here there must be a parent node
        parent_node = self.get_best_parent(col, row)

        score_change = self.get_score_change(col, row, parent_node, match)
        # score_change = 1
        # print("Col: {}, Row: {}. Letters: {}, {}\nMatch: {}, parent at ({}, {}), weight of {}".format(col, row, self.list1[col], self.list2[row], match, parent_node[COL], parent_node[ROW], self.get(parent_node).weight))
        return Node(score_change, parent_node, match)

    def get_best_parent(self, col, row):
        left = self.get_left(col, row)
        up = self.get_up(col, row)
        diag = self.get_diag(col, row)
        direction = ''
        node = ''

        if diag and self.get(diag).weight >= self.get_max(left, up):
            node = diag
            direction = 'diag'
        elif left and not up:
            node = left
            direction = 'left'
        elif up and not left:
            node = up
            direction = 'up'
        else:
            if self.get(left).weight > self.get(up).weight:
                node = left
                direction = 'left'
            else:
                node = up
                direction = 'up'
        return {'node':node, 'direction':direction}

    def get(self, coord):
        if coord is None or coord[0] < 0 or coord[1] < 0:
            return None
        return self.dag_array[coord[COL]][coord[ROW]]

    def get_max(self, first, second):
        return max(self.get(first).weight, self.get(second).weight)

    def get_left(self, col, row):
        if col == 0:
            return None
        return (col - 1, row)
        # return self.dag_array[col - 1][row]

    def get_up(self, col, row):
        if row == 0:
            return None
        return (col, row - 1)
        # return self.dag_array[col][row - 1]

    def get_diag(self, col, row):
        if row == 0 or col == 0:
            return None
        return (col - 1, row - 1)
        # return self.dag_array[col - 1][row - 1]

    def get_score_change(self, col, row, parent_node, match):
        node_val = 0
        if match:
            node_val += self.match_weight
        else: 
            node_val -= self.miss_weight
        if parent_node['direction'] is 'left' or parent_node['direction'] is 'up':
            node_val -= self.indel
        #else:
            #num_misses = self.num_misses(parent_node) + 1
            #score = round(-(match_weight / num_misses))
            #if score < 0:
            #  return 0
            #return score
        node_val += self.get(parent_node['node']).weight
        if node_val < 0:
          return 0
        return node_val

    def num_misses(self, coords):
        node = self.get(coords)
        if node.parent is None:
            return 0

        if not self.get(node.parent).match:
            return self.num_misses(node.parent) + 1
        else:
            return 0

    def get_string(self, coords):
        return ''
        node = self.get(coords)

        if node.match:
            string_array = [self.list1[coords[ROW]]]
        else:
            string_array = [self.NON_MATCH_CHAR]

        string_array = self.get_previous_string(node.parent) + string_array

        return ''.join(string_array)

    def get_previous_string(self, coords):
        node = self.get(coords)

        if not node:
            return None

        if node.match:
            string_array = [self.list1[coords[ROW]]]
        else:
            string_array = [self.NON_MATCH_CHAR]

        previous_string_array = self.get_previous_string(node.parent)

        # If not previous_string_array, we've reached the end
        if previous_string_array:
            string_array = previous_string_array + string_array

        return string_array

    def print_array(self):
        for row in self.dag_array:
            print([node.weight for node in row])

        print("Bottom right: {}".format(self.get_string((self.num_rows - 1, self.num_rows - 1))))

    def print_contigs(self):
        self.print_contigs_from_list(self.list_of_contigs)

    def print_contigs_from_list(self, list_of_contigs):
        print "\nContigs found:"
        for contig in list_of_contigs:
            if contig:
                #print contig
                text = [self.list2[y] for x,y in contig]
                text.reverse()
                print text
        
    def print_score(self):
        score = 0
        for contig in self.list_of_contigs:
            if contig:
                score += self.get(contig[0]).weight
        print "\nScore (total match score): ",# divided by string len): ",
        print score #/ ((len(self.list1) +len(self.list2)) / 2)

    def save_csv(self, filename):
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow([''] + self.list2)
            for index, row in enumerate(self.dag_array):
                writer.writerow([self.list1[index]] + [node.weight for node in row])

    def create_best_non_overlapping_path(self):
        # line up contigs for ref strand (x axis)
        sorted_contigs = sorted(self.list_of_contigs, key = lambda tup: tup[-1][-1])
        contig_dag = [contig_node(contig, self.get(contig[0]).weight) for contig in sorted_contigs]
        # for each contig, 
        for contig in contig_dag:
            # look back and pick the best contig, and add the score to this
            best_parent = self.get_best_parent_contig(contig_dag, contig)
            if best_parent: 
                contig.weight += best_parent.weight
                contig.parent_node = best_parent
        # highest contig score is winner, print it and its path
        best_contig = max(contig_dag, key=lambda cont: cont.weight)
        self.print_contig_path(best_contig)
        return

    def get_best_parent_contig(self, sorted_contigs, this_contig):
        possible_parents = self.get_possible_parent_contigs(sorted_contigs, this_contig)
        best_parent = None
        if possible_parents:
            best_parent = max(possible_parents, key=lambda par: par.weight)
        return best_parent

    def get_possible_parent_contigs(self, sorted_contigs, this_contig):
        # skip any contigs which overlap on cur_strand (ex: 1 "the" matching 2 places)
        possible_parents = [a_contig for a_contig in sorted_contigs if a_contig.contig_path[0][-1] < this_contig.contig_path[-1][-1]]
        #possible_parents = self.remove_overlapping_contigs(possible_parents, this_contig)
        return possible_parents

   # def remove_overlapping_contigs(self, possible_contigs, matching_contig):
   #     for pos_contig in possible_parents:
   #     return

   # def get_first_and_last(self, contig):
   #     return
    #def make_chain_of_contigs_from_class

    def print_contig_path(self, contig):
        print "\nBest contig path score: " + str(contig.weight)
        temp_contig = contig
        contigs_list = []
        if not contig:
            return
        while temp_contig:
            #print temp_contig.contig_path
            contigs_list.append(temp_contig.contig_path)
            temp_contig = temp_contig.parent_node
        self.print_contigs_from_list(contigs_list)
        return

if __name__ == '__main__':
    #string1 = list("The mouse is a house, and I ate it")
    #string2 = list("I ate my louse, and bought a mouse")
    reference = list("""aaaggggcttgatttatttcccatttttcttaaattcctcctttttagagagccaggccactctccttactgagatgcccaa""")
    patient = list("""aagccaactgagatgcccaaagggggccactctccttgcttttcctcctttttagaggatttatttcccatttttcttaaaa""")

    dag = DAG(reference, patient)

    dag.compute_all_long_contigs()

    dag.print_contigs()
    dag.print_score()
    dag.create_best_non_overlapping_path()
    dag.save_csv("dag.csv")
