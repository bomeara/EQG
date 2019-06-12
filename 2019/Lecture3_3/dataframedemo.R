iris			# A built-in dataframe
head(iris)	# just first 6 rows
dim(iris)	# dimensions

iris[1,2]	# accessing value 1,2
iris[1:3, 1]	# pay attention to which return vectors
iris[3, 4:1]	# and which return dataframes, and why
iris[1, "Species"]

iris["Species"]
head(iris["Species"])
head(iris[4:5])

iris$Species

class(iris["Species"])
class(iris$Species)
dim(iris["Species"])
dim(iris$Species)
length(iris$Species)

iris[[5]]
class(iris[[5]])

length(iris)
