from pprint import pprint
import csv

COL = 0
ROW = 1

match_weight = 10
NON_MATCH_CHAR = '#'


class Node(object):
    """docstring for Node"""
    def __init__(self, weight, parent, match):

        self.weight = weight
        # parent is (x,y) coord
        self.parent = parent

        self.match = match


class DAG(object):
    """docstring for DAG"""
    def __init__(self, list1, list2):
        self.list1 = list1
        self.list2 = list2

        self.num_rows = len(self.list1)
        self.num_cols = len(self.list2)

        self.dag_array = None

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
                return Node(match_weight, None, match)
            else:
                return Node(-match_weight, None, match)

        # If we get here there must be a parent node
        parent_node = self.get_best_parent(col, row)

        score_change = self.get_score_change(col, row, parent_node, match)
        # score_change = 1
        # print("Col: {}, Row: {}. Letters: {}, {}\nMatch: {}, parent at ({}, {}), weight of {}".format(col, row, self.list1[col], self.list2[row], match, parent_node[COL], parent_node[ROW], self.get(parent_node).weight))

        return Node(self.get(parent_node).weight + score_change, parent_node, match)

    def get_best_parent(self, col, row):
        left = self.get_left(col, row)
        up = self.get_up(col, row)
        diag = self.get_diag(col, row)

        if diag and self.get(diag).weight >= self.get_max(left, up):
            return diag
        elif left and not up:
            return left
        elif up and not left:
            return up
        else:
            if self.get(left).weight > self.get(up).weight:
                return left
            else:
                return up

    def get(self, coord):
        if coord is None or coord[0] < 0 or coord[1] < 0:
            return None

        return self.dag_array[coord[COL]][coord[ROW]]

    def get_max(self, first, second):
        return max(self.get(first).weight, self.get(second).weight)

        # if first and second:
        #     return max(self.get(first).weight, self.get(second).weight)
        # elif first and not second:
        #     return self.get(first).weight
        # else:
        #     return self.get(second).weight

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

        if match:
            return match_weight
        else:
            num_misses = self.num_misses(parent_node) + 1
            return round(-(match_weight / num_misses))

            # if parent is diagnal
            # if parent_node[COL] < col and parent_node[ROW] < row:
            #     return match_weight
            # else:
            #     return 0

    def num_misses(self, coords):
        node = self.get(coords)
        if node.parent is None:
            return 0

        if not self.get(node.parent).match:
            return self.num_misses(node.parent) + 1
        else:
            return 0

    def get_string(self, coords):
        node = self.get(coords)

        if node.match:
            string_array = [self.list1[coords[ROW]]]
        else:
            string_array = [NON_MATCH_CHAR]

        string_array = self.get_previous_string(node.parent) + string_array

        return ''.join(string_array)

    def get_previous_string(self, coords):
        node = self.get(coords)

        if not node:
            return None

        if node.match:
            string_array = [self.list1[coords[ROW]]]
        else:
            string_array = [NON_MATCH_CHAR]

        previous_string_array = self.get_previous_string(node.parent)

        # If not previous_string_array, we've reached the end
        if previous_string_array:
            string_array = previous_string_array + string_array

        return string_array

    def print_array(self):
        for row in self.dag_array:
            print([node.weight for node in row])

        print("Bottom right: {}".format(self.get_string((self.num_rows - 1, self.num_rows - 1))))

    def save_csv(self, filename):
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow([''] + self.list2)
            for index, row in enumerate(self.dag_array):
                writer.writerow([self.list1[index]] + [node.weight for node in row])

if __name__ == '__main__':
    string1 = list("The mouse is a house, and I ate it")
    string2 = list("I ate my louse, and bought a mouse")
    # string1 = string2
    # string1 = list("Tree ")
    # string2 = list("Three")

    dag = DAG(string1, string2)

    dag.compare()

    dag.print_array()
    dag.save_csv("dag.csv")
