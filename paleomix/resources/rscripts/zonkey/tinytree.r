#!/usr/bin/env Rscript
# Copyright (c) 2015 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Required for 'read.tree'
library(ape)
library(ggplot2)
library(grid)
library(methods)


TTBar <- setRefClass("TTBar",
            fields = list(leftmax = "numeric",
                          left = "numeric",
                          right = "numeric",
                          rightmax = "numeric"))


TTNode <- setRefClass("TTNode",
            fields = list(children = "list",
                          bar = 'TTBar',
                          len = "numeric",
                          label = "character",
                          group = "character"
                          ),
            methods = list(
                "initialize" = function(children = NULL, bar = NULL, len = 0,
                                        label = "", group = "") {
                    .self$children <- as.list(children)
                    if (!is.null(bar)) {
                        .self$bar <- bar
                    }
                    .self$len <- len
                    .self$label <- label
                    .self$group <- group
                },

                "show" = function() {
                    print(node$pformat())
                },

                "pformat" = function() {
                    return(paste(to_str(), ";", sep=""))
                },

                "to_str" = function() {
                    fields <- NULL
                    if (length(children)) {
                        child_str <- NULL
                        for (child in children) {
                            child_str <- c(child_str, child$to_str())
                        }

                        fields <- c(fields, "(", paste(child_str, sep="", collapse=","), ")")
                    }

                    if (nchar(label) > 0) {
                        fields <- c(fields, label)
                    }

                    if (length(len) > 0) {
                        fields <- c(fields, ":", len)
                    }

                    return(paste(fields, sep="", collapse=""))
                },

                "height_above" = function() {
                    total <- ifelse(length(children) > 0, 0, 1)
                    for (child in children_above()) {
                        total <- total + child$height_above() + child$height_below()
                    }

                    return(total)
                },

                "height_below" = function() {
                    total <- ifelse(length(children) > 0, 0, 1)
                    for (child in children_below()) {
                        total <- total + child$height_above() + child$height_below()
                    }

                    return(total)
                },

                "height" = function() {
                    return(height_above() + height_below())
                },

                "width" = function() {
                    total <- 0

                    for (child in children) {
                        total <- max(total, child$width())
                    }

                    if (length(len) > 0) {
                        total <- total + len
                    }

                    return(total)
                },

                "to_tables" = function(from_x=0, from_y=0) {
                    current_x <- from_x + ifelse(length(len) > 0, len, 0)

                    # Horizontal line
                    tables <- list(
                        labels=data.frame(
                            start_x=current_x,
                            start_y=from_y + calc_offset(from_y),
                            label=label,
                            group=get_labelgroup()
                        ),
                        segments=data.frame(
                            start_x=from_x,
                            start_y=from_y + calc_offset(from_y),
                            end_x=current_x,
                            end_y=from_y + calc_offset(from_y),
                            group=get_linegroup()),
                        bars=to_bar(current_x, from_y + calc_offset(from_y)))

                    max_y <- from_y
                    current_y <- max_y
                    for (child in children_above()) {
                        current_y <- current_y + child$height_below()
                        tables <- merge_tables(tables, child$to_tables(current_x, current_y))
                        max_y <- current_y + child$calc_offset(current_y)
                        current_y <- current_y + child$height_above()
                    }

                    min_y <- from_y
                    current_y <- min_y
                    for (child in children_below()) {
                        current_y <- current_y - child$height_above()
                        tables <- merge_tables(tables, child$to_tables(current_x, current_y))
                        min_y <- current_y + child$calc_offset(current_y)
                        current_y <- current_y - child$height_below()
                    }

                    # Vertical line
                    tables$segments <- rbind(tables$segments,
                        data.frame(
                            start_x=current_x,
                            start_y=max_y,
                            end_x=current_x,
                            end_y=min_y,
                            group=get_linegroup()))

                    return(tables)
                },

                "to_bar" = function(current_x, current_y) {
                    if (length(c(bar$leftmax, bar$left, bar$right, bar$rightmax)) != 4) {
                        return(data.frame())
                    }

                    return(data.frame(
                        start_x=c(bar$leftmax, bar$left, bar$right) + current_x,
                        end_x=c(bar$left, bar$right, bar$rightmax) + current_x,
                        start_y=current_y - c(0.25, 0.5, 0.25),
                        end_y=current_y + c(0.25, 0.5, 0.25)))
                },

                "merge_tables" = function(tbl_a, tbl_b) {
                    result <- list()
                    for (name in unique(names(tbl_a), names(tbl_b))) {
                        result[[name]] <- rbind(tbl_a[[name]], tbl_b[[name]])
                    }
                    return(result)
                },

                "clade" = function(taxa) {
                    # FIXME: Handle multiple tips with identical label
                    if (length(intersect(taxa, tips())) != length(taxa)) {
                        return(NULL)
                    }

                    for (child in children) {
                        if (length(intersect(taxa, child$tips())) == length(taxa)) {
                            return(child$clade(taxa))
                        }
                    }

                    return(.self)
                },

                "tips" = function(taxa) {
                    if (length(children) == 0) {
                        return(label)
                    } else {
                        result <- NULL
                        for (child in children) {
                            result <- c(result, child$tips())
                        }
                        return(result)
                    }
                },

                "calc_offset" = function(from_y) {
                    max_y <- from_y
                    current_y <- from_y
                    for (child in children_above()) {
                        current_y <- current_y + child$height_below()
                        max_y <- current_y + child$calc_offset(current_y)
                        current_y <- current_y + child$height_above()
                    }

                    min_y <- from_y
                    current_y <- from_y
                    for (child in children_below()) {
                        current_y <- current_y - child$height_above()
                        min_y <- current_y + child$calc_offset(current_y)
                        current_y <- current_y - child$height_below()
                    }

                    return(max_y - from_y - (max_y - min_y) / 2)
                },

                "children_above" = function() {
                    if (length(children) < 1) {
                        return(list())
                    }
                    return(children[1:ceiling(length(children) / 2)])
                },

                "children_below" = function() {
                    if (length(children) < 1) {
                        return(list())
                    }
                    return(children[(ceiling(length(children) / 2) + 1):length(children)])
                },

                "get_labelgroup" = function(prefix=NULL) {
                    if (is.null(prefix)) {
                        prefix <- ifelse(length(children) > 0, "node", "leaf")
                    }

                    if (nchar(group) > 0) {
                        prefix <- paste(prefix, ":", group, sep="")
                    }

                    return(prefix)
                },

                "get_linegroup" = function() {
                    prefix <- "line"
                    if (nchar(group) > 0) {
                        prefix <- paste(prefix, ":", group, sep="")
                    }

                    return(prefix)
                },

                "set_group" = function(value=NULL) {
                    .self$group <- ifelse(is.null(value), "line", value)
                    for (child in children) {
                        child$set_group(value)
                    }
                },

                "is_leaf" = function() {
                    '
                    Convinience function; returns true if the node is a leaf.
                    '
                    return(length(children) == 0)
                },

                "collect" = function() {
                    '
                    Returns a vector of all in the tree, including this node.
                    '
                    result <- .self
                    for (child in children) {
                        result <- c(result, child$collect())
                    }
                    return(result)
                },

                "sort_nodes" = function() {
                    if (!is_leaf()) {
                        widths <- NULL
                        for (child in children) {
                            widths <- c(widths, child$width())
                            child$sort_nodes()
                        }

                        .self$children <- children[order(widths, decreasing=FALSE)]
                    }
                }))


print.TTNode <- function(node)
{
    print(node$pformat())
}


tinytree.phylo.to.tt <- function(phylo)
{
    nnodes <- nrow(phylo$edge) + 1
    lengths <- phylo$edge.length
    to.node <- phylo$edge[, 2]
    from.node <- phylo$edge[, 1]

    nodes <- list()
    labels <- c(phylo$tip.label, phylo$node.label)
    for (edge in 1:nnodes) {
        len <- lengths[to.node == edge]
        nodes[[edge]] <- TTNode(label = labels[edge],
                                len = as.numeric(len))
    }

    for (edge in 1:nnodes) {
        from <- from.node[to.node == edge]
        if (length(from) != 0 && from != 0) {
            children <- nodes[[from]]$children
            children[[length(children) + 1]] <- nodes[[edge]]
            nodes[[from]]$children <- children
        }
    }

    root <- nodes[[length(phylo$tip.label) + 1]]
    root$len <- 0

    return(root)
}


tinytree.read.newick <- function(filename)
{
	return(tinytree.phylo.to.tt(read.tree(filename)))
}


tinytree.defaults.collect <- function(tt, defaults, values)
{
    stopifnot(!any(is.null(names(values))) || length(values) == 0)

    # Overwrite using user supplied values
    for (idx in seq(values)) {
        defaults[[names(values)[idx]]] <- values[[idx]]
    }

    # Set default values based on type (line, node, leaf, etc.)
    for (node in tt$collect()) {
        for (type in c(node$get_labelgroup(), node$get_linegroup())) {
            if (!(type %in% names(defaults))) {
                root <- unlist(strsplit(type, ":"))[[1]]
                stopifnot(root %in% names(defaults))
                defaults[[type]] <- defaults[[root]]
            }
        }
    }

    return(defaults)
}


tinytree.default.colours <- function(pp, tt, ...)
{
    defaults <- c("line"="black",
                  "node"="darkgrey",
                  "leaf"="black",
                  "bar"="blue")
    defaults <- tinytree.defaults.collect(tt, defaults, list(...))

    return(pp +
           scale_colour_manual(values=defaults) +
           scale_fill_manual(values=defaults))
}


tinytree.default.sizes <- function(pp, tt, ...)
{
    defaults <- c("line"=0.5,
                  "node"=4,
                  "leaf"=5)
    defaults <- tinytree.defaults.collect(tt, defaults, list(...))

    return(pp + scale_size_manual(values=defaults))
}


tinytree.draw <- function(tt, default.scales=TRUE, xaxis="scales", padding=0.3)
{
    tbl <- tt$to_tables(-tt$len)

    pp <- ggplot()
    pp <- pp + geom_segment(data=tbl$segments, lineend="round",
                            aes(x=start_x, y=start_y, xend=end_x, yend=end_y,
                                color=group, size=group))

    if (nrow(tbl$bars) > 0) {
        pp <- pp + geom_rect(data=tbl$bars, alpha=0.3,
                             aes(xmin=start_x, xmax=end_x, ymin=start_y, ymax=end_y,
                                 fill="bar"))
    }

    if (any(!is.na(tbl$labels$label))) {
        labels <- tbl$labels[!is.na(tbl$labels$label),]
        pp <- pp + geom_text(data=labels, hjust=0,
                             aes(label=sprintf(" %s", label),
                                 x=start_x, y=start_y, color=group, size=group))
    }

    pp <- pp + theme_minimal()

    # Disable legend
    pp <- pp + theme(legend.position="none",
    # Disable y axis + y axis labels + grid
                     axis.ticks.y=element_blank(),
                     axis.text.y=element_blank(),
                     panel.grid.minor.y=element_blank(),
                     panel.grid.major.y=element_blank(),
                     panel.grid.major  = element_line(colour = "grey90", size = 0.4),
                     panel.grid.minor  = element_line(colour = "grey90", size = 0.2))

    if (xaxis != "axis") {
        stopifnot(xaxis %in% c("scales", "none"))
        pp <- pp + theme(axis.ticks.x=element_blank(),
                         axis.text.x=element_blank(),
                         panel.grid.minor.x=element_blank(),
                         panel.grid.major.x=element_blank())

        if (xaxis == "scales") {
            y_offset <- min(tbl$segments$start_y, tbl$segments$end_y) - 3
            x_offset <- max(tbl$segments$end_x) * 0.2

            df <- data.frame(x=0, y=y_offset, xend=x_offset, yend=y_offset)
            pp <- pp + geom_segment(data=df,
                                    aes(color="line", size="line",
                                        x=x, xend=xend, y=y, yend=yend))

            df <- data.frame(x=x_offset, y=y_offset, label=paste("", signif(x_offset, 2)))
            pp <- pp + geom_text(data=df, aes(x=x, y=y, hjust=0, label=label,
                                              size="leaf", colour="leaf"))
        }
    }

    # Disable axis labels by default
    pp <- pp + xlab(NULL)
    pp <- pp + ylab(NULL)

    # Default colors; may be overwritten
    if (default.scales) {
        pp <- tinytree.default.sizes(pp, tt)
        pp <- tinytree.default.colours(pp, tt)
    }

    range <- max(tbl$segments$end_x) - min(tbl$segments$start_x)
    pp <- pp + coord_cartesian(xlim=c(min(tbl$segments$start_x),
                                      max(tbl$segments$end_x) + padding * range))

    return(pp)
}


plot.tree <- function(filename, sample_names, padding=0.3)
{
    samples <- read.table(sample_names, as.is=TRUE, comment.char="", header=TRUE)
    tt <- tinytree.read.newick(filename)
    tt$sort_nodes()

    for (node in tt$collect()) {
        if (node$is_leaf()) {
            node$set_group(node$label)
        }
    }

    pp <- tinytree.draw(tt,
                        default.scales=FALSE,
                        padding=padding)
    pp <- tinytree.default.sizes(pp, tt, "node"=3, "leaf"=4, "line"=0.75)

    defaults <- c("line"="black",
                  "node"="grey40",
                  "leaf"="black",
                  "bar"="blue")
    defaults <- tinytree.defaults.collect(tt, defaults, list())

    for (row in 1:nrow(samples)) {
        row <- samples[row, , drop=FALSE]
        key <- sprintf("leaf:%s", row$Name)
        print(c(key, row$Color))

        defaults[[key]] <- row$Color
    }
    print(defaults)

    return(pp +
           scale_colour_manual(values=defaults) +
           scale_fill_manual(values=defaults))
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    cat("Usage: ggtinytree.R <input_file> <sample_names> <output_prefix>\n", file=stderr())
    quit(status=1)
}

input_file <- args[1]
sample_names <- args[2]
output_prefix <- args[3]

pdf(paste(output_prefix, ".pdf", sep=""))
plot.tree(input_file, sample_names)
dev.off()

# bitmap is preferred, since it works in a headless environment
bitmap(paste(output_prefix, ".png", sep=""), height=6, width=6, res=96, taa=4, gaa=4)
plot.tree(input_file, sample_names)
dev.off()

