// phylo.spec.js is created with Cypress, https://www.cypress.io/
// direct download of Cypress is required to run the tests

describe('Parse Newick tree', () => {
    var df = []

    it('readTree: Parse a Newick tree string into a doubly-linked list of JS Objects', () => {
        // Correct return value is dependent on rectLayout() and fortify()
        cy.fixture('../fixtures/timetree.txt').as('timetree')

        df = [{'parentId': 16, 'parentLabel': 'NODE_0000001', 'thisId': 0, 'thisLabel': 'B.1.116', 'children': [], 'branchLength': 0.12477, 'angle': undefined, 'x': 0.12477, 'y': 0, "isTip": true},
        {'parentId': 15, 'parentLabel': 'NODE_0000002', 'thisId': 1, 'thisLabel': 'B.1.422', 'children': [], 'branchLength': 0.06739, 'angle': undefined, 'x': 0.06739, 'y': 1, "isTip": true},
        {'parentId': 4, 'parentLabel': 'NODE_0000003', 'thisId': 2, 'thisLabel': 'B.1.576', 'children': [], 'branchLength': 0, 'angle': undefined, 'x': 0.06193, 'y': 2, "isTip": true},
        {'parentId': 4, 'parentLabel': 'NODE_0000003', 'thisId': 3, 'thisLabel': 'B.1.363', 'children': [], 'branchLength': 0.29781, 'angle': undefined, 'x': 0.35974, 'y': 3, "isTip": true},
        {'parentId': 14, 'parentLabel': 'NODE_0000000', 'thisId': 4, 'thisLabel': 'NODE_0000003', 'children': [2,3], 'branchLength': 0.0364, 'angle': undefined, 'x': 0.06193, 'y': 2.5, "isTip": false},
        {'parentId': 7, 'parentLabel': 'NODE_0000006', 'thisId': 5, 'thisLabel': 'B.1.264', 'children': [], 'branchLength': 0, 'angle': undefined, 'x': 0.06466, 'y': 4, "isTip": true},
        {'parentId': 7, 'parentLabel': 'NODE_0000006', 'thisId': 6, 'thisLabel': 'B.1.264.1', 'children': [], 'branchLength': 0.81975, 'angle': undefined, 'x': 0.8844099999999999, 'y': 5, "isTip": true},
        {'parentId': 14, 'parentLabel': 'NODE_0000000', 'thisId': 7, 'thisLabel': 'NODE_0000006', 'children': [5,6], 'branchLength': 0.03913, 'angle': undefined, 'x': 0.06466, 'y': 4.5, "isTip": false},
        {'parentId': 10, 'parentLabel': 'NODE_0000007', 'thisId': 8, 'thisLabel': 'B.1.513', 'children': [], 'branchLength': 0, 'angle': undefined, 'x': 0.042800000000000005, 'y': 6, "isTip": true},
        {'parentId': 10, 'parentLabel': 'NODE_0000007', 'thisId': 9, 'thisLabel': 'B.1.433', 'children': [], 'branchLength': 0.03825, 'angle': undefined, 'x': 0.08105000000000001, 'y': 7, "isTip": true},
        {'parentId': 14, 'parentLabel': 'NODE_0000000', 'thisId': 10, 'thisLabel': 'NODE_0000007', 'children': [8,9], 'branchLength': 0.01727, 'angle': undefined, 'x': 0.042800000000000005, 'y': 6.5, "isTip": false},
        {'parentId': 13, 'parentLabel': 'NODE_0000008', 'thisId': 11, 'thisLabel': 'B.1.119', 'children': [], 'branchLength': 0.01639, 'angle': undefined, 'x': 0.056459999999999996, 'y': 8, "isTip": true},
        {'parentId': 13, 'parentLabel': 'NODE_0000008', 'thisId': 12, 'thisLabel': 'B.1.371', 'children': [], 'branchLength': 0, 'angle': undefined, 'x': 0.04007, 'y': 9, "isTip": true},
        {'parentId': 14, 'parentLabel': 'NODE_0000000', 'thisId': 13, 'thisLabel': 'NODE_0000008', 'children': [11,12], 'branchLength': 0.01454, 'angle': undefined, 'x': 0.04007, 'y': 8.5, "isTip": false},
        {'parentId': 15, 'parentLabel': 'NODE_0000002', 'thisId': 14, 'thisLabel': 'NODE_0000000', 'children': [4,7,10,13], 'branchLength': 0.02553, 'angle': undefined, 'x': 0.02553, 'y': 5.5, "isTip": false},
        {'parentId': 16, 'parentLabel': 'NODE_0000001', 'thisId': 15, 'thisLabel': 'NODE_0000002', 'children': [1,14], 'branchLength': 0, 'angle': undefined, 'x': 0, 'y': 3.25, "isTip": false},
        {'parentId': null, 'parentLabel': null, 'thisId': 16, 'thisLabel': 'NODE_0000001', 'children': [0,15], 'branchLength': 0, 'angle': undefined, 'x': 0, 'y': 1.625, "isTip": false}]

        cy.fixture('timetree').then((timetree) => {
            cy.window().then((win) => {
                expect(win.readTree(timetree)).to.eql(df)
            })
        })
    })
    it('edges: Generate edge list with x,y coordinates extracted from the respective nodes', () => {
        
        let edges = [
            {'x1': 0.12477, 'y1': 0, 'id1': 0, 'x2': 0, 'y2': 1.625, 'id2': 16, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.06739, 'y1': 1, 'id1': 1, 'x2': 0, 'y2': 3.25, 'id2': 15, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.06193, 'y1': 2, 'id1': 2, 'x2': 0.06193, 'y2': 2.5, 'id2': 4, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.35974, 'y1': 3, 'id1': 3, 'x2': 0.06193, 'y2': 2.5, 'id2': 4, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.06193, 'y1': 2.5, 'id1': 4, 'x2': 0.02553, 'y2': 5.5, 'id2': 14, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.06466, 'y1': 4, 'id1': 5, 'x2': 0.06466, 'y2': 4.5, 'id2': 7, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.8844099999999999, 'y1': 5, 'id1': 6, 'x2': 0.06466, 'y2': 4.5, 'id2': 7, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.06466, 'y1': 4.5, 'id1': 7, 'x2': 0.02553, 'y2': 5.5, 'id2': 14, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.042800000000000005, 'y1': 6, 'id1': 8, 'x2': 0.042800000000000005, 'y2': 6.5, 'id2': 10, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.08105000000000001, 'y1': 7, 'id1': 9, 'x2': 0.042800000000000005, 'y2': 6.5, 'id2': 10, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.042800000000000005, 'y1': 6.5, 'id1': 10, 'x2': 0.02553, 'y2': 5.5, 'id2': 14, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.056459999999999996, 'y1': 8, 'id1': 11, 'x2': 0.04007, 'y2': 8.5, 'id2': 13, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.04007, 'y1': 9, 'id1': 12, 'x2': 0.04007, 'y2': 8.5, 'id2': 13, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.04007, 'y1': 8.5, 'id1': 13, 'x2': 0.02553, 'y2': 5.5, 'id2': 14, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0.02553, 'y1': 5.5, 'id1': 14, 'x2': 0, 'y2': 3.25, 'id2': 15, 'last_date': undefined, 'first_date': undefined},
            {'x1': 0, 'y1': 3.25, 'id1': 15, 'x2': 0, 'y2': 1.625, 'id2': 16, 'last_date': undefined, 'first_date': undefined}
        ]

        cy.window().then((win) => {
            expect(win.edges(df)).to.eql(edges)
        })

    })
    it('traverse: Recursive function for traversal of tree', () => {

    })
})

