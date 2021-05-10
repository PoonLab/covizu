/**
 * Parse a Newick tree string into a doubly-linked
 * list of JS Objects.  Assigns node labels, branch
 * lengths and node IDs (numbering terminal before
 * internal nodes).
 * @param {string} text Newick tree string.
 * @return {object} Root of tree.
 */
function readTree(text) {
    // remove whitespace
    text = text.replace(/ |\t|\r?\n|\r/g, '');

    var tokens = text.split(/(;|\(|\)|,)/),
        root = {'parent': null, 'children':[]},
        curnode = root,
        nodeId = 0,
        nodeinfo;

    var node_labels = [];

    for (const token of tokens) {
        if (token == "" || token == ';') {
            continue
        }
        //console.log(token);
        if (token == '(') {
            // add a child to current node
            var child = {
                'parent': curnode,
                'children': []
            };
            curnode.children.push(child);
            curnode = child;  // climb up
        }
        else if (token == ',') {
            // climb down, add another child to parent
            curnode = curnode.parent;
            var child = {
                'parent': curnode,
                'children': []
            }
            curnode.children.push(child);
            curnode = child;  // climb up
        }
        else if (token == ')') {
            // climb down twice
            curnode = curnode.parent;
            if (curnode === null) {
                break;
            }
        }
        else {
            nodeinfo = token.split(':');
            node_labels.push(nodeinfo);

            if (nodeinfo.length==1) {
                if (token.startsWith(':')) {
                    curnode.label = "";
                    curnode.branchLength = parseFloat(nodeinfo[0]);
                } else {
                    curnode.label = nodeinfo[0];
                    curnode.branchLength = null;
                }
            }
            else if (nodeinfo.length==2) {
                curnode.label = nodeinfo[0];
                curnode.branchLength = parseFloat(nodeinfo[1]);
            }
            else {
                // TODO: handle edge cases with >1 ":"
                console.warn(token, "I don't know what to do with two colons!");
            }
            curnode.id = nodeId++;  // assign then increment
        }
    }

    // if root node is unlabelled
    if (node_labels.length < nodeId) {
        curnode.id = nodeId;
    }

    return (getTimeTreeData(root));
}


/**
 * Recursive function for traversal of tree
 * (output parent before children).
 * @param {object} node
 * @param {string} order: 'preorder' or 'postorder' traversal
 * @param {Array} list: an Array of nodes
 * @return An Array of nodes in pre-order
 */
function traverse(node, order='preorder', list=Array()) {
    if (order=='preorder') list.push(node);
    for (var i=0; i < node.children.length; i++) {
        list = traverse(node.children[i], order, list);
    }
    if (order=='postorder') list.push(node);
    return(list);
}


/**
 * Convert parsed Newick tree from readTree() into more convenient
 * tabular data frame.
 * @param {object} tree: Return value of readTree
 * @param {boolean} sort: if true, sort data frame by node name
 * @return Array of Objects
 */
function fortify(tree, sort=true) {
    var df = [];

    for (const node of traverse(tree, 'preorder')) {
        if (node.parent === null) {
            df.push({
                'parentId': null,
                'parentLabel': null,
                'thisId': node.id,
                'thisLabel': node.label,
                'children': node.children.map(x=>x.id),
                'branchLength': 0.,
                'isTip': (node.children.length===0),
                'x': node.x,
                'y': node.y,
                'angle': node.angle
            })
        }
        else {
            df.push({
                'parentId': node.parent.id,
                'parentLabel': node.parent.label,
                'thisId': node.id,
                'thisLabel': node.label,
                'children': node.children.map(x=>x.id),
                'branchLength': node.branchLength,
                'isTip': (node.children.length===0),
                'x': node.x,
                'y': node.y,
                'angle': node.angle
            })
        }
    }

    if (sort) {
        df = df.sort(function(a, b) {
            return a.thisId - b.thisId;
        })
    }
    return(df);
}


/**
 * Generate edge list with x,y coordinates extracted from the
 * respective nodes.
 * @param {Array} df:  tabular data frame from fortify()
 * @param {boolean} rectangular:  if true, draw two line segments connected
 *                                by right angle
 * @returns {Array}
 */
function edges(df, rectangular=false) {
    var result = [],
        parent, pair;

    // make sure data frame is sorted
    df.sort(function(a, b) {
        return a.thisId - b.thisId;
    });

    for (const row of df) {
        if (row.parentId === null) {
            continue  // skip the root
        }
        parent = df[row.parentId];
        if (parent === null || parent === undefined) {
            console.log('parent null/undefined');
            continue;
        }

        if (rectangular) {
          pair = {
              x1: row.x, y1: row.y, id1: row.thisId,
              x2: parent.x, y2: row.y, id2: undefined,
              last_date: row.isTip ? row.last_date : undefined,
              first_date: row.isTip ? row.first_date : undefined

          };
          result.push(pair);
          pair = {
              x1: parent.x, y1: row.y, id1: undefined,
              x2: parent.x, y2: parent.y, id2: row.parentId,
              last_date: row.isTip ? row.last_date : undefined,
              first_date: row.isTip ? row.first_date : undefined

          };
          result.push(pair);
        }
        else {
          pair = {
              x1: row.x, y1: row.y, id1: row.thisId,
              x2: parent.x, y2: parent.y, id2: row.parentId,
              last_date: row.isTip ? row.last_date : undefined,
              first_date: row.isTip ? row.first_date : undefined
          };
          result.push(pair);
        }
    }
    return(result);
}

