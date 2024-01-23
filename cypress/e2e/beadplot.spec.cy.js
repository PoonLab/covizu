import clusters from '../../data_test/clusters.json'

describe('Beadplot components', () => {
    it('Edge Labels', ()=> {
        cy.window().then((win) => {
            win.reset_tree(true);
            const tip_id = win.display_id.non_recombinants.last - 2;
            cy.get(`[id=id-${tip_id}]`).trigger('click');
            cy.get('.beadplot-content>svg>g>text').should(($text) => {
                expect($text).to.have.length(14)
                expect($text.eq(0)).to.contain('USA/AZ-CDC-LC0932884')
                expect($text.eq(1)).to.contain('USA/AZ-CDC-STM-TNFPTQHEE')
                expect($text.eq(2)).to.contain('USA/AZ-CDC-LC0951048')
            }); 
        });
    })
    it('Samples', () => {
        cy.get('.beadplot-content>svg>g>circle').then(($circ) => {
            expect($circ).to.have.length(23)
            for (let i = 0; i < 8; i++) {
                cy.get($circ.eq(i)).should('have.attr', 'cy', 20)
            }
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
        var tips = {'allregions': {'North America': 23}, 'country': {'USA': 23}}

        cy.window().then((win) => {
            win.gentable(tips)
        })
        cy.get('table:visible>tbody>tr').should('have.length', 1)
        cy.get('table:visible>tbody>tr').eq(0).contains('North America')
        cy.get('table:visible>tbody>tr').eq(0).contains('USA')
        cy.get('table:visible>tbody>tr').eq(0).contains('23')
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

        let targetValue = 5.01
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
        cy.get('[stroke="#bbd"]').should('have.length', 12)
    })
})

describe('Tooltips', () => {
    it('Horizontal Edge', () => {
        cy.get('#USA-AZ-CDC-LC0956766').first().trigger('mouseover')
        cy.get('.tooltip').contains('Parent: USA/AZ-CDC-LC0932884')
        cy.get('.tooltip').contains('North America: 2')
        cy.get('.tooltip').contains('Unique collection dates: 2')
        cy.get('.tooltip').contains('2022-12-05 / 2022-12-12')
    })
    it('Bead', () => {
        cy.get('#EPI_ISL_1048746').first().trigger('mouseover')
        cy.get('.tooltip').contains('Parent: USA/AZ-CDC-LC0932884')
        cy.get('.tooltip').contains('Genomic distance: 4.06')
        cy.get('.tooltip').contains('North America: 1')
        cy.get('.tooltip').contains('Collection date: 2022-12-05')
    })
})