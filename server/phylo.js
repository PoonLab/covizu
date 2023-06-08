const {$DATA_FOLDER} = require("../config")
const dbstats = require(`../${$DATA_FOLDER}/dbstats.json`)

/**
 * Parse a Newick tree string into a doubly-linked
 * list of JS Objects.  Assigns node labels, branch
 * lengths and node IDs (numbering terminal before
 * internal nodes).
 * @param {string} text Newick tree string.
 * @return {object} Root of tree.
 */
const readTree = (text) => {
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
 * Get the data frame
 * @param {Object} timetree:  time-scaled phylogenetic tree imported as JSON
 */
const getTimeTreeData = (timetree) => {
  // generate tree layout (x, y coordinates
  rectLayout(timetree);

  var df = fortify(timetree);

  return(df);
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
        'rawLabel': dbstats["lineages"][node.label] ? dbstats["lineages"][node.label]["raw_lineage"] : null,
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
      'rawLabel': dbstats["lineages"][node.label] ? dbstats["lineages"][node.label]["raw_lineage"] : null,
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
 * Rectangular layout of tree, update nodes in place with x,y coordinates
 * @param {object} root
 */
function rectLayout(root) {
  // assign vertical positions to tips by postorder traversal
  var counter = 0;
  for (const node of traverse(root, 'postorder')) {
    if (node.children.length === 0) {
      // assign position to tip
      node.y = counter;
      counter++;
    } else {
      // ancestral node position is average of child nodes
      node.y = 0;
      for (var i = 0; i < node.children.length; i++) {
        var child = node.children[i];
        node.y += child.y;
      }
      node.y /= node.children.length;
    }
  }

  // assign horizontal positions by preorder traversal
  for (const node of traverse(root, 'preorder')) {
    if (node.parent === null) {
      // assign root to x=0
      node.x = 0.;
    } else {
      node.x = node.parent.x + node.branchLength;
    }
  }
}


module.exports = {
  readTree
};