// covizu.spec.js is created with Cypress, https://www.cypress.io/
// direct download of Cypress is required to run the tests

describe('CoVizu page', () => {
    it('Sucessfully loads', () => {
        cy.visit("/")
    })
    it('Language options redirect', () => {

        cy.contains('a', 'es').click()
        cy.url().should('include', '/index-es.html')

        cy.contains('a', 'fr').click()
        cy.url().should('include', '/index-fr.html')

        cy.contains('a', 'zh').click()
        cy.url().should('include', '/index-zh.html')

        cy.contains('a', 'en').click()
        cy.url().should('include', '/index.html')
    })
    it('Welcome dialog opens automatically', () => {
        cy.get('#splash').should('be.visible')
        cy.get('#splash-button').should('be.disabled')
    })
    it('Splash button is enabled after clusters.json loads', () => {
        cy.get('#splash-button', {timeout: 40000}).should('be.enabled')
        cy.get('#splash-button').click().should('not.be.visible')
    })
})

describe('Help dialogs', () => {
    it('Open and close', () => {
        cy.get('[style="cursor: help"]').each(($el) => {
            cy.get($el).click()
            cy.get('.ui-dialog:visible').find('.ui-dialog-titlebar-close').click()
            cy.get('.ui-dialog').should('not.be.visible')
        })
    })
})

describe('Search Interface', () => {
    it('Buttons are initially disabled', () => {
        cy.get('#navigation>button').each(($el) => {
            cy.get($el).should('be.disabled')
        })
    })
    it('Error message for invalid date formats', () => {
        cy.get('#start-date').type('2021-05-13')
        cy.get('#end-date').type('2021-05-10')
        cy.get('#search-button').click({force:true})
        cy.contains('Start Date must be before the End Date.')

        cy.get('#start-date').clear()
        cy.get('#end-date').clear()
        cy.get('#start-date').type('20210510')
        cy.get('#search-button').click({force:true})
        cy.contains('Invalid Start Date (YYYY-MM-DD)')

        cy.get('#start-date').clear()
        cy.get('#end-date').clear()
        cy.get('#end-date').type('20210513')
        cy.get('#search-button').click({force:true})
        cy.contains('Invalid End Date (YYYY-MM-DD)')
        cy.get('#end-date').clear()
    })
    it("Searching 'HMH' results in selection and expected point count", () => {
        cy.get('#search-input').type('HMH')
        cy.get('#search-button').click({force:true})
        cy.get('.selectionH').should('exist')
        cy.get('.SelectedCluster').should('exist')
        cy.window().then((win) => {
            cy.get('#search_stats').contains(`1 of ${win.bead_id_to_accession.length} points`)
        })
    })
    it("Clicking the next button highlights the correct bead", () => {
        cy.window().then((win) => {
            cy.get('.selectionH').should('have.attr', 'bead', `${win.bead_id_to_accession[0]}`)
            cy.get('#next_button').click()
            cy.get('.selectionH').should('have.length', 1)
            cy.get('.selectionH').should('have.attr', 'bead', `${win.bead_id_to_accession[1]}`)
        })
    })
    it('Clear button resets search interface and cluster selection', () => {
        cy.get('#clear_button').click()
        cy.get('#search-bar>input').each(($el) => {
            cy.get($el).should('be.empty')
        })
        cy.get('#search_stats').contains('0 of 0 points')
        cy.get('.selectionH').should('not.exist')
        cy.get('.SelectedCluster').should('not.exist')
    })
})

describe('Tooltips', () => {
    it('Appear on hover over cluster', () => {
        cy.get('rect:visible').first().trigger('mouseover')
        cy.get('.tooltip').should('be.visible').should('have.length', 1)
    })
    it('Cluster tooltips contain relevant and correct information', () => {
        cy.get('rect:visible').first().trigger('mouseover').invoke('attr','id').as('id')
        cy.get('@id').then(id => {
            var id_number = id.substring(3)
            cy.window().then((win) => {
                cy.get('.tooltip').contains(`Number of variants: ${win.tips[id_number]['varcount']}`)
                cy.wrap(win.tips[id_number]['mutations']).each(($el) => {
                    cy.get('.tooltip').contains($el)
                })
                cy.wrap(Object.keys(win.tips[id_number]['allregions'])).each(($el, index) => {
                    cy.get('.tooltip').contains(`${$el}: ${Object.values(win.tips[id_number]['allregions'])[index]}`)
                })        
            })
        })
    })
    it('Appear on hover over bead', () => {
        cy.get('circle:visible').first().trigger('mouseover')
        cy.get('.tooltip').should('be.visible').should('have.length', 1)
    })
})