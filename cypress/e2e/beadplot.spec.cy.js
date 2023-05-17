import clusters from '../../data_test/clusters.json'

describe('Beadplot components', () => {
    it('Edge Labels', ()=> {
        cy.get('.beadplot-content>svg>g>text').should(($text) => {
            expect($text).to.have.length(3)
            expect($text.eq(0)).to.contain('unsampled0')
            expect($text.eq(1)).to.contain('Scotland/GCVR-16F930')
            expect($text.eq(2)).to.contain('Scotland/GCVR-170556')
        }) 
    })
    it('Samples', () => {
        cy.get('.beadplot-content>svg>g>circle').then(($circ) => {
            expect($circ).to.have.length(4)
            for (let i = 0; i < 3; i++) {
                cy.get($circ.eq(i)).should('have.attr', 'cy', 45)
            }
            cy.get($circ.eq(3)).should('have.attr', 'cy', 70)
        }) 

    })
})

describe('Tables', () => {
    it('region_to_string', () => {
        cy.window().then((win) => {
            expect(win.region_to_string({'North America':5,'Asia':10})).to.eq("<b>Number of cases:</b><br>&nbsp;&nbsp;North America: 5<br>&nbsp;&nbsp;Asia: 10<br>Total: 15<br>")
        })
    })
    it('Country details: gentable', () => {
        cy.get('#tabs').contains('Countries').click()
        var tips = {'allregions': {'North America': 13}, 'country': {'Canada': 8, 'USA': 5}}

        cy.window().then((win) => {
            win.gentable(tips)
        })
        cy.get('table:visible>tbody>tr').should('have.length', 2)
        cy.get('table:visible>tbody>tr').eq(0).contains('North America')
        cy.get('table:visible>tbody>tr').eq(0).contains('Canada')
        cy.get('table:visible>tbody>tr').eq(0).contains('8')

        cy.get('table:visible>tbody>tr').eq(1).contains('North America')
        cy.get('table:visible>tbody>tr').eq(1).contains('USA')
        cy.get('table:visible>tbody>tr').eq(1).contains('5')
    })
})

describe('Edge slider', () => {
    var targetValue, currentValue, steps, arrows = 0
    const increment = 0.01

    it('Arrows increment slider by 0.01', () => {
        cy.get('#left-arrow').click()
        cy.get('#custom-handle').contains('1.99')

        cy.get('#right-arrow').click()
        cy.get('#custom-handle').contains('2')
    })
    it('Max slider value is equal to max edge length in cluster', () => {

        // https://stackoverflow.com/questions/64855669/moving-slider-with-cypress

        let targetValue = 12
        let currentValue = 2
        let increment = 0.01
        steps = (targetValue - currentValue) / increment + 1
        arrows = '{rightarrow}'.repeat(steps)

        cy.get('#custom-handle').contains(currentValue).type(arrows)
        cy.get('#custom-handle', { timeout: 7000 }).should('contain', targetValue)

        // Number displayed should not increase once max is reached
        cy.get('#custom-handle').type('{rightarrow}')
        cy.get('#custom-handle').should('contain', targetValue)
    })
    it('All edges are visible on max slider value ', () => {
        cy.get('[stroke="#bbd"]').should('have.length', 2)
    })
    it('No edges are visible on slider value of 11.99 ', () => {
        cy.get('#left-arrow').click()
        cy.get('#custom-handle').should('contain', 11.99)
        cy.get('[stroke="#bbd"]').should('not.exist')
    })
})

describe('Tooltips', () => {
    it('Horizontal Edge', () => {
        cy.get('#Scotland-GCVR-16F930').first().trigger('mouseover')
        cy.get('.tooltip').contains('Parent: unsampled0')
        cy.get('.tooltip').contains('Europe: 3')
        cy.get('.tooltip').contains('Unique collection dates: 3')
        cy.get('.tooltip').contains('2020-03-19 / 2020-03-27')
    })
    it('Bead', () => {
        cy.get('circle:visible').first().trigger('mouseover')
        cy.get('.tooltip').contains('Parent: unsampled0')
        cy.get('.tooltip').contains('Genomic distance: 12')
        cy.get('.tooltip').contains('Europe: 1')
        cy.get('.tooltip').contains('Collection date: 2020-03-19')
    })
})