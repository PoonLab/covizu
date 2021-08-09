// drawtree.spec.js is created with Cypress, https://www.cypress.io/
// direct download of Cypress is required to run the tests

import clusters from '../fixtures/clusters.json'
import countries from '../fixtures/countries.json'
import dbstats from '../fixtures/dbstats.json'

describe('Selected Lineage', () => {
    it('Default open rectangle around first lineage', () => {
        cy.get('rect').last().should('have.attr', 'class', 'clicked')
    })
    it('Draw an open rectangle around a clicked lineage', () => {
        cy.get('rect:visible').first().click().invoke('attr','id').as('id')
        cy.get('@id').then(id => {
            cy.get(`#${id}`).should('have.attr', 'class', 'clicked')
            cy.get('rect.clicked').should('have.length', 1)
            cy.get('text.clicked').should('have.length', 1)
        })
    })
})

describe('drawtree', () => {
    var df = []

    it('No negative branch lengths', () => {
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
    })

    cy.window().then((win) => {
        win.drawtree(df)
    })
})