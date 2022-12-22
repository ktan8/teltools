class InputListLine():
        '''
        Class to reprsent what the input
        list line looks like
        '''
        def __init__(self, line):
                self.line = line
                self.lineArr = line.strip().split("\t")
                self.chrom = self.lineArr[0]
                self.position = self.lineArr[1]
                self.dist_left = None
                self.dist_right = None
                self.near_edge = None

        def update_dist_variables(self, dist_left, dist_right):
                '''
                Update the variables for the
                distance from left and right
                '''
                self.dist_left = dist_left
                self.dist_right = dist_right

        def update_nearedge(self, near_edge):
                '''
                Update whether the site is near
                the edge or not.
                '''
                self.near_edge = near_edge

        def get_output(self):
                '''
                Generate and return the output as a list
                '''
                output = self.lineArr.append(self.dist_left, self.dist_right, self.near_edge)

                return output

class InputList():
        '''
        Input list of sites that can be extracted
        and analyzed
        '''
        def __init__(self, sitefile):
                self.sitefile = sitefile

        def calc_edge_distance(self, fai_dict):
                '''
                Calculate how far the
                '''
                for line_class in self.site_file_generator():
                        dist_left_edge = line_class.position
                        dist_right_edge = fai_dict[line_class.chrom] - line_class.position

                        line_class.update_dist(dist_left_edge, dist_right_edge)

                        yield line_class


        def site_file_generator(self):
                '''
                Generator for the site file
                '''
                f = open(self.sitefile, "r")

                for line in f:
                        line_class = inputListLine
                        yield line_class


        def check_near_edge(self, edge_cutoff, fai_dict):
                '''
                Check if the site is near the edge
                or not
                '''
                for line_class in self.calc_edge_distance(fai_dict):
                        near_edge = None
                        if line_class.dist_left < edge_cutoff or \
                        line_class.dist_right < edge_cutoff:
                                near_edge = 1
                        else:
                                near_edge = 0


class SiteFile:
	def __init__(self, site_file):
		self.site_file = site_file
		

	def line_generator(self, outputList = False):
		'''
		Returns each of the line in the file
		'''
		f = open(self.site_file, "r")
		
		for line in f:
			if outputList:
				line_list = line.strip().split("\t")
				yield line_list
			else:
				yield line.strip()
		
	

