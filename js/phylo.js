
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
    text = text.replace(/ \t/g, '');

    var tokens = text.split(/(;|\(|\)|,)/),
        root = {'parent': null, 'children':[]},
        curnode = root,
        nodeId = 0;

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
            var nodeinfo = token.split(':');

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
    curnode.id = nodeId;

    return (root);
    //return({json: root, error: null});
}

//var s = "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;";
//readTree(s);

/**
 * Recursive function for pre-order traversal of tree
 * (output parent before children).
 * @param {object} node
 * @param {Array} list An Array of nodes
 * @return An Array of nodes in pre-order
 */
function preorder(node, list=[]) {
    list.push(node);
    for (var i=0; i < node.children.length; i++) {
        list = preorder(node.children[i], list);
    }
    return(list);
}

function postorder(node, list=[]) {
    for (var i=0; i < node.children.length; i++) {
        list = postorder(node.children[i], list);
    }
    list.push(node);
    return(list);
}

function levelorder(root) {
    // aka breadth-first search
    var queue = [root],
        result = [],
        curnode;

    while (queue.length > 0) {
        curnode = queue.pop();
        result.push(curnode);
        for (const child of curnode.children) {
            queue.push(child);
        }
    }
    return(result);
}

/**
 * Count the number of tips that descend from this node
 * @param {object} thisnode
 */
function numTips(thisnode) {
    var result = 0;
    for (const node of levelorder(thisnode)) {
        if (node.children.length == 0) result++;
    }
    return(result);
}


/**
 * Convert parsed Newick tree from readTree() into data
 * frame.
 * @param {object} tree Return value of readTree
 * @return Array of Objects
 */
function fortify(tree, sort=true) {
    var df = [];

    for (const node of preorder(tree)) {
        if (node.parent === null) {
            df.push({
                'parentId': null,
                'parentLabel': null,
                'thisId': node.id,
                'thisLabel': node.label,
                'children': node.children.map(x=>x.id),
                'branchLength': 0.,
                'isTip': (node.children.length==0),
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
                'isTip': (node.children.length==0),
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


function edges(df, rectangular=false) {
    var result = [],
        parent, pair;

    // make sure data frame is sorted
    df.sort(function(a, b) {
        return a.thisId - b.thisId;
    })

    for (const row of df) {
        x1 = row.x;
        y1 = row.y;
        if (row.parentId === null) {
            continue  // skip the root
        }
        parent = df[row.parentId];
        if (parent === null || parent === undefined) continue;

        if (rectangular) {
          var pair = {
              x1: row.x, y1: row.y, id1: row.thisId,
              x2: parent.x, y2: row.y, id2: undefined
          };
          result.push(pair);
          var pair = {
              x1: parent.x, y1: row.y, id1: undefined,
              x2: parent.x, y2: parent.y, id2: row.parentId
          };
          result.push(pair);
        }
        else {
          var pair = {
              x1: row.x, y1: row.y, id1: row.thisId,
              x2: parent.x, y2: parent.y, id2: row.parentId
          };
          result.push(pair);
        }
    }
    return(result);
}


/**
 * Equal-angle layout algorithm for unrooted trees.
 * Populates the nodes of a tree object with information on
 * the angles to draw branches such that they do not
 * intersect.
 * @param {object} node
 */
function equalAngleLayout(node) {
    if (node.parent === null) {
        // node is root
        node.start = 0.;  // guarantees no arcs overlap 0
        node.end = 2.; // *pi
        node.angle = 0.;  // irrelevant
        node.ntips = numTips(node);
        node.x = 0;
        node.y = 0;
    }

    var child, arc, lastStart = node.start;

    for (var i=0; i<node.children.length; i++) {
        child = node.children[i];
        child.ntips = numTips(child);

        // assign proportion of arc to this child
        arc = (node.end-node.start) * child.ntips/node.ntips;
        child.start = lastStart;
        child.end = child.start + arc;

        // bisect the arc
        child.angle = child.start + (child.end-child.start)/2.;
        lastStart = child.end;

        // map to coordinates
        child.x = node.x + child.branchLength * Math.sin(child.angle*Math.PI);
        child.y = node.y + child.branchLength * Math.cos(child.angle*Math.PI);

        // climb up
        equalAngleLayout(child);
    }
}


/**
 * Rectangular layout
 * @param {object} root
 */
function rectLayout(root) {

}
