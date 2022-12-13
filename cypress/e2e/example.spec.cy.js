/**
these tests were copied from the covizu.spec.cy.js but the failing tests were 
removed to check the accuracy of the Github Actions CI/CD.
*/

describe('CoVizu page', () => {
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
})
